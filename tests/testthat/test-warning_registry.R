# test-warning_registry.R
# Comprehensive unit tests for lwdid warning registry mechanism
# Covers: constructor, register, flush, verbose levels, aggregation,
#         get_log, get_diagnostics, count/clear/merge, end-to-end, post-flush

# ============================================================================
# Group 1: Constructor tests
# ============================================================================

test_that("new_warning_registry() returns an environment", {
  reg <- new_warning_registry()
  expect_true(is.environment(reg))
})

test_that("initial records is an empty list", {
  reg <- new_warning_registry()
  expect_identical(reg$records, list())
})

test_that("initial flushed is FALSE", {
  reg <- new_warning_registry()
  expect_false(reg$flushed)
})

test_that("constructor returns environment with emptyenv() parent", {
  reg <- new_warning_registry()
  expect_identical(parent.env(reg), emptyenv())
})

test_that("registry contains all 7 method functions", {
  reg <- new_warning_registry()
  methods <- c("register", "flush", "get_log",
               "get_diagnostics", "count", "clear", "merge")
  for (m in methods) {
    expect_true(is.function(reg[[m]]),
                info = paste("Method", m, "should be a function"))
  }
})

# ============================================================================
# Group 2: register method tests
# ============================================================================

test_that("register buffers records without producing warnings", {
  reg <- new_warning_registry()
  expect_silent(
    reg$register("lwdid_small_sample", "Sample size below 30",
                 cohort = 2005L, period = 2003L)
  )
  expect_equal(length(reg$records), 1L)
})

test_that("register accumulates multiple records", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg1")
  reg$register("lwdid_convergence", "msg2")
  reg$register("lwdid_small_sample", "msg3")
  expect_equal(length(reg$records), 3L)
})

test_that("register rejects non-character category (numeric)", {
  reg <- new_warning_registry()
  expect_error(
    reg$register(123, "msg"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("register rejects multi-element category vector", {
  reg <- new_warning_registry()
  expect_error(
    reg$register(c("a", "b"), "msg"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("register rejects non-character message (NULL)", {
  reg <- new_warning_registry()
  expect_error(
    reg$register("lwdid_small_sample", NULL),
    class = "lwdid_invalid_parameter"
  )
})

test_that("register rejects non-character message (numeric)", {
  reg <- new_warning_registry()
  expect_error(
    reg$register("lwdid_small_sample", 42),
    class = "lwdid_invalid_parameter"
  )
})

test_that("context=NULL is converted to list() internally", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg", context = NULL)
  expect_identical(reg$records[[1]]$context, list())
  # get_diagnostics should work without error
  diag <- reg$get_diagnostics()
  expect_true(is.list(diag))
})

test_that("each record contains a numeric timestamp", {
  reg <- new_warning_registry()
  before <- as.numeric(Sys.time())
  reg$register("lwdid_small_sample", "msg")
  after <- as.numeric(Sys.time())
  ts <- reg$records[[1]]$timestamp
  expect_true(is.numeric(ts))
  expect_true(ts >= before && ts <= after)
})

test_that("register() returns invisible(NULL)", {
  reg <- new_warning_registry()
  expect_invisible(reg$register("lwdid_small_sample", "msg"))
  result <- reg$register("lwdid_small_sample", "msg2")
  expect_null(result)
})

# ============================================================================
# Group 3: flush method tests
# ============================================================================

test_that("flush is idempotent: second flush produces no warnings", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg1")
  # First flush should produce warnings
  expect_warning(reg$flush("default"), "msg1")
  # Second flush should be silent
  expect_silent(reg$flush("default"))
  # Records still accessible via get_log
  log_df <- reg$get_log()
  expect_equal(nrow(log_df), 1L)
})

test_that("flush rejects invalid verbose value", {
  reg <- new_warning_registry()
  expect_error(
    reg$flush(verbose = "invalid"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("flush verbose is case-insensitive", {
  # Test "Default"
  reg1 <- new_warning_registry()
  reg1$register("lwdid_small_sample", "msg")
  expect_warning(reg1$flush(verbose = "Default"))

  # Test "QUIET"
  reg2 <- new_warning_registry()
  reg2$register("lwdid_convergence", "msg")
  expect_warning(reg2$flush(verbose = "QUIET"))

  # Test "Verbose"
  reg3 <- new_warning_registry()
  reg3$register("lwdid_small_sample", "msg")
  expect_warning(reg3$flush(verbose = "Verbose"))
})

test_that("flush validates verbose BEFORE checking flushed flag", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg")
  reg$flush("default")  # flush once, now flushed=TRUE
  # Even though already flushed, invalid verbose should still error
  expect_error(
    reg$flush(verbose = "invalid"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("flush on empty registry: no output, flushed becomes TRUE", {
  reg <- new_warning_registry()
  expect_silent(reg$flush("default"))
  expect_true(reg$flushed)
})

# ============================================================================
# Group 4: verbose three-level control tests
# ============================================================================

test_that("quiet mode: only critical categories are emitted", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "small sample msg")
  reg$register("lwdid_convergence", "convergence msg")

  # Collect all warnings
  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("quiet"),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # Only lwdid_convergence (critical) should be emitted
  expect_equal(length(warnings_caught), 1L)
  expect_true(inherits(warnings_caught[[1]], "lwdid_convergence"))
  # lwdid_small_sample should NOT appear
  classes <- vapply(warnings_caught, function(w) class(w)[1],
                    character(1))
  expect_false("lwdid_small_sample" %in% classes)
})

test_that("default mode: all categories emitted as aggregated summaries", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "small msg",
               cohort = 2005L, period = 2003L)
  reg$register("lwdid_small_sample", "small msg",
               cohort = 2006L, period = 2004L)
  reg$register("lwdid_convergence", "conv msg",
               cohort = 2005L, period = 2003L)
  reg$register("lwdid_convergence", "conv msg",
               cohort = 2007L, period = 2005L)

  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default"),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # 2 categories -> 2 aggregated summaries
  expect_equal(length(warnings_caught), 2L)
  # Each should have [category] prefix format
  msgs <- vapply(warnings_caught, conditionMessage, character(1))
  expect_true(all(grepl("^\\[", msgs)))
})

test_that("verbose mode: each record emitted individually", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg A")
  reg$register("lwdid_convergence", "msg B")
  reg$register("lwdid_small_sample", "msg C")

  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("verbose"),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # 3 records -> 3 individual warnings
  expect_equal(length(warnings_caught), 3L)
  msgs <- vapply(warnings_caught, conditionMessage, character(1))
  expect_equal(msgs, c("msg A", "msg B", "msg C"))
  # Verbose mode: no [category] prefix
  expect_false(any(grepl("^\\[", msgs)))
})

test_that("flush output preserves category for tryCatch", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "small sample issue",
               cohort = 2005L, period = 2003L)

  caught <- NULL
  withCallingHandlers(
    reg$flush("default"),
    lwdid_small_sample = function(w) {
      caught <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_false(is.null(caught))
  expect_true(inherits(caught, "lwdid_small_sample"))
  expect_true(inherits(caught, "lwdid_warning"))
})

# ============================================================================
# Group 5: aggregation logic tests
# ============================================================================

test_that("aggregation: M/T ratio format with total_pairs", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "small sample",
               cohort = 2005L, period = 2003L)
  reg$register("lwdid_small_sample", "small sample",
               cohort = 2006L, period = 2004L)

  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default", total_pairs = 120L),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  msg <- conditionMessage(warnings_caught[[1]])
  # Should contain "2/120 (cohort, period) pairs"
  expect_true(grepl("2/120 (cohort, period) pairs", msg,
                    fixed = TRUE))
})

test_that("aggregation: N pairs format without total_pairs", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "small sample",
               cohort = 2005L, period = 2003L)
  reg$register("lwdid_small_sample", "small sample",
               cohort = 2006L, period = 2004L)

  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default", total_pairs = NULL),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  msg <- conditionMessage(warnings_caught[[1]])
  # Should contain "2 (cohort, period) pairs" without ratio
  expect_true(grepl("2 (cohort, period) pairs", msg,
                    fixed = TRUE))
  expect_false(grepl("/", msg))
})

test_that("aggregation: N occurrences format without cohort/period", {
  reg <- new_warning_registry()
  reg$register("lwdid_numerical", "numerical issue")
  reg$register("lwdid_numerical", "numerical issue")
  reg$register("lwdid_numerical", "numerical issue")

  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default"),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  msg <- conditionMessage(warnings_caught[[1]])
  expect_true(grepl("3 occurrences", msg, fixed = TRUE))
})

test_that("aggregation: total_pairs=0 treated as NULL", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg",
               cohort = 2005L, period = 2003L)

  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default", total_pairs = 0L),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  msg <- conditionMessage(warnings_caught[[1]])
  # Should use N pairs format, not M/T ratio
  expect_true(grepl("1 (cohort, period) pairs", msg,
                    fixed = TRUE))
  expect_false(grepl("/0", msg))
})

test_that("aggregation: non-integer total_pairs (120.0) works", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg",
               cohort = 2005L, period = 2003L)

  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default", total_pairs = 120.0),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  msg <- conditionMessage(warnings_caught[[1]])
  expect_true(grepl("1/120 (cohort, period) pairs", msg,
                    fixed = TRUE))
})

test_that("aggregation: total_pairs>0 but no cohort/period -> occurrences", {
  reg <- new_warning_registry()
  reg$register("lwdid_numerical", "numerical issue")
  reg$register("lwdid_numerical", "numerical issue")

  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default", total_pairs = 120L),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  msg <- conditionMessage(warnings_caught[[1]])
  # n_affected=0, so should use occurrences format
  expect_true(grepl("2 occurrences", msg, fixed = TRUE))
  expect_false(grepl("/120", msg))
})

# ============================================================================
# Group 6: get_log method tests
# ============================================================================

test_that("get_log returns data.frame with 5 columns", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg",
               cohort = 2005L, period = 2003L)
  log_df <- reg$get_log()
  expect_s3_class(log_df, "data.frame")
  expect_equal(ncol(log_df), 5L)
  expect_equal(names(log_df),
               c("category", "message", "cohort", "period", "timestamp"))
})

test_that("get_log works before and after flush with same result", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg",
               cohort = 2005L, period = 2003L)
  log_before <- reg$get_log()
  suppressWarnings(reg$flush("default"))
  log_after <- reg$get_log()
  expect_equal(log_before, log_after)
})

test_that("get_log returns zero-row data.frame for empty registry", {
  reg <- new_warning_registry()
  log_df <- reg$get_log()
  expect_s3_class(log_df, "data.frame")
  expect_equal(nrow(log_df), 0L)
  expect_equal(ncol(log_df), 5L)
  expect_equal(names(log_df),
               c("category", "message", "cohort", "period", "timestamp"))
})

test_that("get_log: NULL cohort/period become NA_integer_", {
  reg <- new_warning_registry()
  reg$register("lwdid_numerical", "msg")  # no cohort/period
  log_df <- reg$get_log()
  expect_true(is.na(log_df$cohort[1]))
  expect_true(is.integer(log_df$cohort))
  expect_true(is.na(log_df$period[1]))
  expect_true(is.integer(log_df$period))
})

test_that("get_log: mixed NULL and non-NULL cohort/period", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg1",
               cohort = 2005L, period = 2003L)
  reg$register("lwdid_numerical", "msg2")  # NULL cohort/period
  log_df <- reg$get_log()
  expect_equal(nrow(log_df), 2L)
  expect_equal(log_df$cohort[1], 2005L)
  expect_true(is.na(log_df$cohort[2]))
  expect_equal(log_df$period[1], 2003L)
  expect_true(is.na(log_df$period[2]))
})

test_that("get_log: timestamp column is numeric", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg")
  log_df <- reg$get_log()
  expect_true(is.numeric(log_df$timestamp))
})

# ============================================================================
# Group 7: get_diagnostics method tests
# ============================================================================

test_that("get_diagnostics returns aggregated summary with 5 fields", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg1",
               cohort = 2005L, period = 2003L,
               context = list(n_treated = 3))
  reg$register("lwdid_small_sample", "msg1",
               cohort = 2006L, period = 2004L,
               context = list(n_treated = 5))
  reg$register("lwdid_convergence", "msg2",
               cohort = 2005L, period = 2003L,
               context = list(iterations = 25))
  diag <- reg$get_diagnostics()
  expect_equal(length(diag), 2L)
  # Each entry has 5 fields
  for (entry in diag) {
    expect_true(all(c("category", "message", "count",
                      "affected_pairs", "context_summary")
                    %in% names(entry)))
  }
})

test_that("get_diagnostics returns empty list for empty registry", {
  reg <- new_warning_registry()
  expect_identical(reg$get_diagnostics(), list())
})

test_that("get_diagnostics: affected_pairs only includes non-NULL", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg",
               cohort = 2005L, period = 2003L)
  reg$register("lwdid_small_sample", "msg")  # no cohort/period
  reg$register("lwdid_small_sample", "msg",
               cohort = 2006L, period = 2004L)
  diag <- reg$get_diagnostics()
  # 3 records, but only 2 have cohort/period
  expect_equal(length(diag[[1]]$affected_pairs), 2L)
})

test_that("get_diagnostics: NULL->NA_integer_ in affected_pairs", {
  reg <- new_warning_registry()
  # cohort non-NULL but period NULL
  reg$register("lwdid_small_sample", "msg",
               cohort = 2005L, period = NULL)
  diag <- reg$get_diagnostics()
  pair <- diag[[1]]$affected_pairs[[1]]
  expect_equal(length(pair), 2L)
  expect_equal(pair[1], 2005L)
  expect_true(is.na(pair[2]))
  expect_true(is.integer(pair[2]))
})

test_that("get_diagnostics: context_summary keys sorted alphabetically", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg",
               context = list(z_score = 1.5, alpha = 0.05,
                              n_treated = 10))
  reg$register("lwdid_small_sample", "msg",
               context = list(z_score = 2.0, alpha = 0.10,
                              n_treated = 20))
  diag <- reg$get_diagnostics()
  ctx <- diag[[1]]$context_summary
  ctx_names <- names(ctx)
  # alpha < n_treated < z_score (alphabetical)
  expected_names <- c("alpha_min", "alpha_max",
                      "n_treated_min", "n_treated_max",
                      "z_score_min", "z_score_max")
  expect_equal(ctx_names, expected_names)
  # Verify values
  expect_equal(ctx$alpha_min, 0.05)
  expect_equal(ctx$alpha_max, 0.10)
  expect_equal(ctx$n_treated_min, 10)
  expect_equal(ctx$n_treated_max, 20)
  expect_equal(ctx$z_score_min, 1.5)
  expect_equal(ctx$z_score_max, 2.0)
})

test_that("get_diagnostics: non-numeric context fields skipped", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg",
               context = list(method = "logit", n_treated = 10))
  reg$register("lwdid_small_sample", "msg",
               context = list(method = "probit", n_treated = 20))
  diag <- reg$get_diagnostics()
  ctx <- diag[[1]]$context_summary
  # "method" should not appear
  expect_false("method_min" %in% names(ctx))
  expect_false("method_max" %in% names(ctx))
  # n_treated should appear
  expect_equal(ctx$n_treated_min, 10)
  expect_equal(ctx$n_treated_max, 20)
})

test_that("get_diagnostics works before and after flush", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg",
               cohort = 2005L, period = 2003L)
  diag_before <- reg$get_diagnostics()
  suppressWarnings(reg$flush("default"))
  diag_after <- reg$get_diagnostics()
  expect_equal(length(diag_before), length(diag_after))
  expect_equal(diag_before[[1]]$count, diag_after[[1]]$count)
})

# ============================================================================
# Group 8: count/clear/merge method tests
# ============================================================================

test_that("count: total and per-category", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg1")
  reg$register("lwdid_small_sample", "msg2")
  reg$register("lwdid_convergence", "msg3")
  expect_equal(reg$count(), 3L)
  expect_equal(reg$count("lwdid_small_sample"), 2L)
  expect_equal(reg$count("lwdid_convergence"), 1L)
})

test_that("count: empty registry returns 0L", {
  reg <- new_warning_registry()
  result <- reg$count()
  expect_equal(result, 0L)
  expect_true(is.integer(result))
})

test_that("count: nonexistent category returns 0L", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg")
  result <- reg$count("nonexistent_category")
  expect_equal(result, 0L)
  expect_true(is.integer(result))
})

test_that("clear: resets records and flushed flag", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg")
  suppressWarnings(reg$flush("default"))
  expect_true(reg$flushed)
  expect_equal(length(reg$records), 1L)

  reg$clear()
  expect_identical(reg$records, list())
  expect_false(reg$flushed)
})

test_that("clear then register and flush works", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg1")
  suppressWarnings(reg$flush("default"))
  reg$clear()
  reg$register("lwdid_convergence", "msg2")
  # Should be able to flush again
  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default"),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(length(warnings_caught), 1L)
  expect_true(inherits(warnings_caught[[1]], "lwdid_convergence"))
})

test_that("merge: combines records from two registries", {
  reg1 <- new_warning_registry()
  reg2 <- new_warning_registry()
  reg1$register("lwdid_small_sample", "msg1")
  reg1$register("lwdid_small_sample", "msg2")
  reg2$register("lwdid_convergence", "msg3")
  reg1$merge(reg2)
  expect_equal(length(reg1$records), 3L)
  expect_equal(reg1$count(), 3L)
})

test_that("merge: resets flushed flag", {
  reg1 <- new_warning_registry()
  reg1$register("lwdid_small_sample", "msg")
  suppressWarnings(reg1$flush("default"))
  expect_true(reg1$flushed)

  reg2 <- new_warning_registry()
  reg2$register("lwdid_convergence", "msg2")
  reg1$merge(reg2)
  expect_false(reg1$flushed)
})

test_that("merge: empty registry is legal, still resets flushed", {
  reg1 <- new_warning_registry()
  reg1$register("lwdid_small_sample", "msg")
  suppressWarnings(reg1$flush("default"))
  expect_true(reg1$flushed)

  reg2 <- new_warning_registry()  # empty
  reg1$merge(reg2)
  expect_equal(length(reg1$records), 1L)  # unchanged
  expect_false(reg1$flushed)  # reset
})

test_that("merge: other_registry is not modified", {
  reg1 <- new_warning_registry()
  reg2 <- new_warning_registry()
  reg2$register("lwdid_convergence", "msg")
  reg2$register("lwdid_numerical", "msg2")
  n_before <- length(reg2$records)

  reg1$merge(reg2)
  expect_equal(length(reg2$records), n_before)
  expect_equal(reg2$records[[1]]$category, "lwdid_convergence")
  expect_equal(reg2$records[[2]]$category, "lwdid_numerical")
})

test_that("merge: non-registry list throws error", {
  reg <- new_warning_registry()
  expect_error(
    reg$merge(list(a = 1)),
    class = "lwdid_invalid_parameter"
  )
})

test_that("merge: non-environment throws error", {
  reg <- new_warning_registry()
  expect_error(
    reg$merge(123),
    class = "lwdid_invalid_parameter"
  )
})

# ============================================================================
# Group 9: end-to-end test (parallel worker simulation)
# ============================================================================

test_that("end-to-end: multi-worker register -> merge -> flush", {
  main_reg <- new_warning_registry()
  worker1 <- new_warning_registry()
  worker2 <- new_warning_registry()

  # Worker 1: small sample warnings
  worker1$register("lwdid_small_sample", "Small sample size",
                   cohort = 2005L, period = 2003L,
                   context = list(n_treated = 3, n_control = 12))
  worker1$register("lwdid_small_sample", "Small sample size",
                   cohort = 2006L, period = 2004L,
                   context = list(n_treated = 5, n_control = 20))

  # Worker 2: convergence + numerical warnings
  worker2$register("lwdid_convergence", "Failed to converge",
                   cohort = 2007L, period = 2005L,
                   context = list(iterations = 100))
  worker2$register("lwdid_numerical", "Singular matrix",
                   context = list(condition_number = 1e15))

  # Merge workers into main
  main_reg$merge(worker1)
  main_reg$merge(worker2)

  # Verify total count
  expect_equal(main_reg$count(), 4L)
  expect_equal(main_reg$count("lwdid_small_sample"), 2L)
  expect_equal(main_reg$count("lwdid_convergence"), 1L)
  expect_equal(main_reg$count("lwdid_numerical"), 1L)

  # Flush with default mode
  warnings_caught <- list()
  withCallingHandlers(
    main_reg$flush("default", total_pairs = 50L),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # 3 categories -> 3 aggregated summaries
  expect_equal(length(warnings_caught), 3L)

  # Verify get_log contains all 4 records
  log_df <- main_reg$get_log()
  expect_equal(nrow(log_df), 4L)
  expect_equal(sum(log_df$category == "lwdid_small_sample"), 2L)
  expect_equal(sum(log_df$category == "lwdid_convergence"), 1L)
  expect_equal(sum(log_df$category == "lwdid_numerical"), 1L)

  # Verify get_diagnostics
  diag <- main_reg$get_diagnostics()
  expect_equal(length(diag), 3L)

  # Find small_sample diagnostics
  ss_diag <- Filter(
    function(d) d$category == "lwdid_small_sample", diag
  )[[1]]
  expect_equal(ss_diag$count, 2L)
  expect_equal(length(ss_diag$affected_pairs), 2L)
  # Context summary: n_control and n_treated
  expect_equal(ss_diag$context_summary$n_control_min, 12)
  expect_equal(ss_diag$context_summary$n_control_max, 20)
  expect_equal(ss_diag$context_summary$n_treated_min, 3)
  expect_equal(ss_diag$context_summary$n_treated_max, 5)

  # Find convergence diagnostics
  conv_diag <- Filter(
    function(d) d$category == "lwdid_convergence", diag
  )[[1]]
  expect_equal(conv_diag$count, 1L)
  expect_equal(conv_diag$context_summary$iterations_min, 100)
  expect_equal(conv_diag$context_summary$iterations_max, 100)

  # Find numerical diagnostics
  num_diag <- Filter(
    function(d) d$category == "lwdid_numerical", diag
  )[[1]]
  expect_equal(num_diag$count, 1L)
  # No cohort/period -> no affected pairs
  expect_equal(length(num_diag$affected_pairs), 0L)
  expect_equal(num_diag$context_summary$condition_number_min, 1e15)
  expect_equal(num_diag$context_summary$condition_number_max, 1e15)
})

# ============================================================================
# Group 10: post-flush behavior tests
# ============================================================================

test_that("post-flush register: new records appended, second flush silent", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "msg1")
  suppressWarnings(reg$flush("default"))

  # Register after flush
  reg$register("lwdid_convergence", "msg2")
  log_df <- reg$get_log()
  expect_equal(nrow(log_df), 2L)

  # Second flush should be silent (idempotent)
  expect_silent(reg$flush("default"))
})

test_that("flush returns invisible(NULL) in all code paths", {
  # Path 1: normal flush with warnings
  reg1 <- new_warning_registry()
  reg1$register("lwdid_small_sample", "msg")
  result1 <- suppressWarnings(reg1$flush("default"))
  expect_null(result1)
  expect_invisible(suppressWarnings({
    reg_tmp <- new_warning_registry()
    reg_tmp$register("lwdid_small_sample", "msg")
    reg_tmp$flush("default")
  }))

  # Path 2: empty registry
  reg2 <- new_warning_registry()
  expect_invisible(reg2$flush("default"))
  result2 <- reg2$flush("default")  # already flushed path
  expect_null(result2)

  # Path 3: already flushed
  reg3 <- new_warning_registry()
  reg3$register("lwdid_small_sample", "msg")
  suppressWarnings(reg3$flush("default"))
  expect_invisible(reg3$flush("default"))
  result3 <- reg3$flush("default")
  expect_null(result3)
})

test_that("register returns invisible(NULL)", {
  reg <- new_warning_registry()
  expect_invisible(reg$register("lwdid_small_sample", "msg"))
  result <- reg$register("lwdid_small_sample", "msg2")
  expect_null(result)
})

test_that("merge flushed registry: target can re-flush", {
  reg_a <- new_warning_registry()
  reg_a$register("lwdid_small_sample", "msg from A",
                 cohort = 2005L, period = 2003L)
  reg_a$register("lwdid_convergence", "conv from A",
                 cohort = 2006L, period = 2004L)
  suppressWarnings(reg_a$flush("default"))
  expect_true(reg_a$flushed)

  reg_b <- new_warning_registry()
  reg_b$merge(reg_a)
  expect_false(reg_b$flushed)

  # B should be able to flush A's records
  warnings_caught <- list()
  withCallingHandlers(
    reg_b$flush("default"),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(length(warnings_caught), 2L)
  log_df <- reg_b$get_log()
  expect_equal(nrow(log_df), 2L)
})

test_that("clear then full workflow: clear->register->flush->get_log", {
  reg <- new_warning_registry()
  reg$register("lwdid_small_sample", "old msg")
  suppressWarnings(reg$flush("default"))

  reg$clear()
  expect_equal(reg$count(), 0L)
  expect_false(reg$flushed)

  reg$register("lwdid_convergence", "new msg",
               cohort = 2010L, period = 2008L)
  warnings_caught <- list()
  withCallingHandlers(
    reg$flush("default"),
    lwdid_warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(length(warnings_caught), 1L)
  expect_true(inherits(warnings_caught[[1]], "lwdid_convergence"))

  log_df <- reg$get_log()
  expect_equal(nrow(log_df), 1L)
  expect_equal(log_df$category[1], "lwdid_convergence")
})

test_that("count return type is always integer", {
  reg <- new_warning_registry()
  # Empty registry
  r1 <- reg$count()
  expect_true(is.integer(r1))

  # With records
  reg$register("lwdid_small_sample", "msg")
  r2 <- reg$count()
  expect_true(is.integer(r2))

  # Per-category
  r3 <- reg$count("lwdid_small_sample")
  expect_true(is.integer(r3))

  # Nonexistent category
  r4 <- reg$count("nonexistent")
  expect_true(is.integer(r4))
})
