# test-conditions.R
# 完整单元测试：lwdid自定义错误与警告条件体系
# 覆盖：工厂函数、便捷抛出函数、S3类链、tryCatch捕获、附加信息字段、消息格式化

# ============================================================================
# 第1组：工厂函数基础测试
# ============================================================================

test_that("lwdid_error() creates correct class chain", {
  cond <- lwdid_error("test msg", "lwdid_test")
  expect_s3_class(cond, "error")
  expect_s3_class(cond, "condition")
  expect_s3_class(cond, "lwdid_error")
  expect_equal(class(cond), c("lwdid_test", "lwdid_error", "error", "condition"))
})

test_that("lwdid_error() supports vector class (intermediate parent)", {
  cond <- lwdid_error("test", c("child", "parent"))
  expect_equal(class(cond),
               c("child", "parent", "lwdid_error", "error", "condition"))
})

test_that("lwdid_error() deduplicates class chain", {
  cond <- lwdid_error("test", c("a", "lwdid_error"))
  expect_equal(class(cond), c("a", "lwdid_error", "error", "condition"))
  # lwdid_error should not appear twice
  expect_equal(sum(class(cond) == "lwdid_error"), 1L)
})

test_that("lwdid_error() stores additional fields", {
  cond <- lwdid_error("msg", "lwdid_test", foo = 42, bar = "baz")
  expect_equal(cond$foo, 42)
  expect_equal(cond$bar, "baz")
})

test_that("lwdid_error() stores call field", {
  my_call <- quote(some_function(x = 1))
  cond <- lwdid_error("msg", "lwdid_test", call = my_call)
  expect_equal(cond$call, my_call)
})

test_that("lwdid_warning() creates correct class chain", {
  cond <- lwdid_warning("test msg", "lwdid_test_warn")
  expect_s3_class(cond, "warning")
  expect_s3_class(cond, "condition")
  expect_s3_class(cond, "lwdid_warning")
  expect_equal(class(cond),
               c("lwdid_test_warn", "lwdid_warning", "warning", "condition"))
})

test_that("lwdid_warning() deduplicates class chain", {
  cond <- lwdid_warning("test", c("a", "lwdid_warning"))
  expect_equal(class(cond), c("a", "lwdid_warning", "warning", "condition"))
  expect_equal(sum(class(cond) == "lwdid_warning"), 1L)
})

test_that("lwdid_warning() stores additional fields", {
  cond <- lwdid_warning("msg", "lwdid_test_warn", x = 10)
  expect_equal(cond$x, 10)
})

# ============================================================================
# 第2组：便捷抛出函数测试
# ============================================================================

test_that("stop_lwdid() throws error condition", {
  expect_error(
    stop_lwdid("boom", "lwdid_test_err"),
    "boom"
  )
})

test_that("stop_lwdid() can be caught by lwdid_error", {
  cond <- tryCatch(
    stop_lwdid("boom", "lwdid_test_err"),
    lwdid_error = function(c) c
  )
  expect_s3_class(cond, "lwdid_error")
  expect_equal(conditionMessage(cond), "boom")
})

test_that("warn_lwdid() issues warning condition", {
  expect_warning(
    warn_lwdid("careful", "lwdid_test_warn"),
    "careful"
  )
})

test_that("warn_lwdid() can be caught by lwdid_warning", {
  cond <- tryCatch(
    {
      warn_lwdid("careful", "lwdid_test_warn")
      NULL
    },
    lwdid_warning = function(c) c
  )
  expect_s3_class(cond, "lwdid_warning")
  expect_equal(conditionMessage(cond), "careful")
})

test_that("stop_lwdid() call field captures caller info", {
  my_func <- function() stop_lwdid("err", "lwdid_test")
  cond <- tryCatch(my_func(), lwdid_error = function(c) c)
  expect_false(is.null(cond$call))
})

test_that("warn_lwdid() call field captures caller info", {
  my_func <- function() warn_lwdid("warn", "lwdid_test")
  cond <- tryCatch(my_func(), lwdid_warning = function(c) c)
  expect_false(is.null(cond$call))
})

# ============================================================================
# 第3组：18种错误条件S3类链测试
# ============================================================================

test_that("lwdid_invalid_parameter class chain is correct", {
  cond <- lwdid_invalid_parameter_error("p", "v", c("a", "b"))
  expect_equal(class(cond),
               c("lwdid_invalid_parameter", "lwdid_error", "error", "condition"))
})

test_that("lwdid_invalid_rolling class chain is correct", {
  cond <- lwdid_invalid_rolling_error("bad")
  expect_equal(class(cond),
               c("lwdid_invalid_rolling", "lwdid_invalid_parameter",
                 "lwdid_error", "error", "condition"))
})

test_that("lwdid_invalid_vce class chain is correct", {
  cond <- lwdid_invalid_vce_error("bad")
  expect_equal(class(cond),
               c("lwdid_invalid_vce", "lwdid_invalid_parameter",
                 "lwdid_error", "error", "condition"))
})

test_that("lwdid_insufficient_data class chain is correct", {
  cond <- lwdid_insufficient_data_error(10, 3, 7)
  expect_equal(class(cond),
               c("lwdid_insufficient_data", "lwdid_error", "error", "condition"))
})

test_that("lwdid_no_treated class chain is correct", {
  cond <- lwdid_no_treated_error("d", 100)
  expect_equal(class(cond),
               c("lwdid_no_treated", "lwdid_insufficient_data",
                 "lwdid_error", "error", "condition"))
})

test_that("lwdid_no_control class chain is correct", {
  cond <- lwdid_no_control_error("d", 100)
  expect_equal(class(cond),
               c("lwdid_no_control", "lwdid_insufficient_data",
                 "lwdid_error", "error", "condition"))
})

test_that("lwdid_insufficient_pre_periods class chain is correct", {
  cond <- lwdid_insufficient_pre_periods_error(1, 2, "demean")
  expect_equal(class(cond),
               c("lwdid_insufficient_pre_periods", "lwdid_insufficient_data",
                 "lwdid_error", "error", "condition"))
})

test_that("lwdid_insufficient_quarter_diversity class chain is correct", {
  cond <- lwdid_insufficient_quarter_diversity_error(1, c(3), c(1, 2), c(3, 4))
  expect_equal(class(cond),
               c("lwdid_insufficient_quarter_diversity",
                 "lwdid_insufficient_data",
                 "lwdid_error", "error", "condition"))
})

test_that("lwdid_no_never_treated class chain is correct", {
  cond <- lwdid_no_never_treated_error("cohort", "never_treated")
  expect_equal(class(cond),
               c("lwdid_no_never_treated", "lwdid_insufficient_data",
                 "lwdid_error", "error", "condition"))
})

test_that("lwdid_unbalanced_panel class chain is correct", {
  cond <- lwdid_unbalanced_panel_error(5, 10, 3, 30.0)
  expect_equal(class(cond),
               c("lwdid_unbalanced_panel", "lwdid_error", "error", "condition"))
})

test_that("lwdid_time_discontinuity class chain is correct", {
  cond <- lwdid_time_discontinuity_error(c(2005, 2008), "state", "year")
  expect_equal(class(cond),
               c("lwdid_time_discontinuity", "lwdid_error", "error", "condition"))
})

test_that("lwdid_missing_column class chain is correct", {
  cond <- lwdid_missing_column_error("y", c("x", "z"))
  expect_equal(class(cond),
               c("lwdid_missing_column", "lwdid_error", "error", "condition"))
})

test_that("lwdid_randomization_error class chain is correct", {
  cond <- lwdid_randomization_error_cond(500, "timeout")
  expect_equal(class(cond),
               c("lwdid_randomization_error", "lwdid_error", "error", "condition"))
})

test_that("lwdid_visualization_error class chain is correct", {
  cond <- lwdid_visualization_error_cond("event_study", "missing data")
  expect_equal(class(cond),
               c("lwdid_visualization_error", "lwdid_error", "error", "condition"))
})

test_that("lwdid_invalid_staggered_data class chain is correct", {
  cond <- lwdid_invalid_staggered_data_error("gvar", c(-1, -2), "negative values")
  expect_equal(class(cond),
               c("lwdid_invalid_staggered_data", "lwdid_error", "error", "condition"))
})

test_that("lwdid_aggregation_error class chain is correct", {
  cond <- lwdid_aggregation_error_cond("overall", "no data")
  expect_equal(class(cond),
               c("lwdid_aggregation_error", "lwdid_error", "error", "condition"))
})

test_that("lwdid_invalid_aggregation class chain is correct", {
  cond <- lwdid_invalid_aggregation_error("event", "violation", "detail")
  expect_equal(class(cond),
               c("lwdid_invalid_aggregation", "lwdid_aggregation_error",
                 "lwdid_error", "error", "condition"))
})

test_that("lwdid_insufficient_cell_size class chain is correct", {
  cond <- lwdid_insufficient_cell_size_error(10, 5, 20)
  expect_equal(class(cond),
               c("lwdid_insufficient_cell_size", "lwdid_aggregation_error",
                 "lwdid_error", "error", "condition"))
})

# ============================================================================
# 第4组：5种警告条件S3类链测试
# ============================================================================

test_that("lwdid_small_sample class chain is correct", {
  cond <- lwdid_small_sample_warning(20, 8, 12)
  expect_equal(class(cond),
               c("lwdid_small_sample", "lwdid_warning", "warning", "condition"))
})

test_that("lwdid_overlap class chain is correct", {
  cond <- lwdid_overlap_warning(5, c(0.01, 0.99), 50.0)
  expect_equal(class(cond),
               c("lwdid_overlap", "lwdid_warning", "warning", "condition"))
})

test_that("lwdid_numerical class chain is correct", {
  cond <- lwdid_numerical_warning("singular matrix", "solve_ols")
  expect_equal(class(cond),
               c("lwdid_numerical", "lwdid_warning", "warning", "condition"))
})

test_that("lwdid_data class chain is correct", {
  cond <- lwdid_data_warning("missing values", "rows dropped")
  expect_equal(class(cond),
               c("lwdid_data", "lwdid_warning", "warning", "condition"))
})

test_that("lwdid_convergence class chain is correct", {
  cond <- lwdid_convergence_warning("logistic", 25, 100, 1e-8)
  expect_equal(class(cond),
               c("lwdid_convergence", "lwdid_warning", "warning", "condition"))
})

# ============================================================================
# 第5组：tryCatch捕获测试
# ============================================================================

test_that("precise catch: specific error condition", {
  cond <- tryCatch(
    stop(lwdid_missing_column_error("y", c("x", "z"))),
    lwdid_missing_column = function(c) c
  )
  expect_s3_class(cond, "lwdid_missing_column")
})

test_that("precise catch: specific warning condition", {
  cond <- tryCatch(
    warning(lwdid_overlap_warning(5, c(0.01, 0.99), 50.0)),
    lwdid_overlap = function(c) c
  )
  expect_s3_class(cond, "lwdid_overlap")
})

test_that("intermediate parent lwdid_insufficient_data catches lwdid_no_treated", {
  cond <- tryCatch(
    stop(lwdid_no_treated_error("d", 100)),
    lwdid_insufficient_data = function(c) c
  )
  expect_s3_class(cond, "lwdid_no_treated")
  expect_s3_class(cond, "lwdid_insufficient_data")
})

test_that("intermediate parent lwdid_insufficient_data catches lwdid_no_control", {
  cond <- tryCatch(
    stop(lwdid_no_control_error("d", 100)),
    lwdid_insufficient_data = function(c) c
  )
  expect_s3_class(cond, "lwdid_no_control")
})

test_that("intermediate parent lwdid_insufficient_data catches lwdid_insufficient_pre_periods", {
  cond <- tryCatch(
    stop(lwdid_insufficient_pre_periods_error(1, 2, "demean")),
    lwdid_insufficient_data = function(c) c
  )
  expect_s3_class(cond, "lwdid_insufficient_pre_periods")
})

test_that("intermediate parent lwdid_insufficient_data catches lwdid_insufficient_quarter_diversity", {
  cond <- tryCatch(
    stop(lwdid_insufficient_quarter_diversity_error(1, c(3), c(1, 2), c(3, 4))),
    lwdid_insufficient_data = function(c) c
  )
  expect_s3_class(cond, "lwdid_insufficient_quarter_diversity")
})

test_that("intermediate parent lwdid_insufficient_data catches lwdid_no_never_treated", {
  cond <- tryCatch(
    stop(lwdid_no_never_treated_error("cohort", "never_treated")),
    lwdid_insufficient_data = function(c) c
  )
  expect_s3_class(cond, "lwdid_no_never_treated")
})

test_that("intermediate parent lwdid_invalid_parameter catches lwdid_invalid_rolling", {
  cond <- tryCatch(
    stop(lwdid_invalid_rolling_error("bad")),
    lwdid_invalid_parameter = function(c) c
  )
  expect_s3_class(cond, "lwdid_invalid_rolling")
})

test_that("intermediate parent lwdid_invalid_parameter catches lwdid_invalid_vce", {
  cond <- tryCatch(
    stop(lwdid_invalid_vce_error("bad")),
    lwdid_invalid_parameter = function(c) c
  )
  expect_s3_class(cond, "lwdid_invalid_vce")
})

test_that("intermediate parent lwdid_aggregation_error catches lwdid_invalid_aggregation", {
  cond <- tryCatch(
    stop(lwdid_invalid_aggregation_error("event", "v", "d")),
    lwdid_aggregation_error = function(c) c
  )
  expect_s3_class(cond, "lwdid_invalid_aggregation")
})

test_that("intermediate parent lwdid_aggregation_error catches lwdid_insufficient_cell_size", {
  cond <- tryCatch(
    stop(lwdid_insufficient_cell_size_error(10, 5, 20)),
    lwdid_aggregation_error = function(c) c
  )
  expect_s3_class(cond, "lwdid_insufficient_cell_size")
})

test_that("top-level lwdid_error catches all 18 error conditions", {
  # Build a list of all 18 error condition constructors with sample args
  error_conds <- list(
    lwdid_invalid_parameter_error("p", "v", c("a")),
    lwdid_invalid_rolling_error("bad"),
    lwdid_invalid_vce_error("bad"),
    lwdid_insufficient_data_error(10, 3, 7),
    lwdid_no_treated_error("d", 100),
    lwdid_no_control_error("d", 100),
    lwdid_insufficient_pre_periods_error(1, 2, "demean"),
    lwdid_insufficient_quarter_diversity_error(1, c(3), c(1, 2), c(3, 4)),
    lwdid_no_never_treated_error("cohort", "never_treated"),
    lwdid_unbalanced_panel_error(5, 10, 3, 30.0),
    lwdid_time_discontinuity_error(c(2005), "state", "year"),
    lwdid_missing_column_error("y", c("x")),
    lwdid_randomization_error_cond(500, "timeout"),
    lwdid_visualization_error_cond("es", "err"),
    lwdid_invalid_staggered_data_error("gvar", c(-1), "neg"),
    lwdid_aggregation_error_cond("overall", "no data"),
    lwdid_invalid_aggregation_error("event", "v", "d"),
    lwdid_insufficient_cell_size_error(10, 5, 20)
  )
  for (i in seq_along(error_conds)) {
    caught <- tryCatch(
      stop(error_conds[[i]]),
      lwdid_error = function(c) c
    )
    expect_true(inherits(caught, "lwdid_error"),
                info = paste("Error condition", i, "not caught by lwdid_error"))
  }
})

test_that("top-level lwdid_warning catches all 5 warning conditions", {
  warn_conds <- list(
    lwdid_small_sample_warning(20, 8, 12),
    lwdid_overlap_warning(5, c(0.01, 0.99), 50.0),
    lwdid_numerical_warning("issue", "func"),
    lwdid_data_warning("issue", "action"),
    lwdid_convergence_warning("logistic", 25, 100, 1e-8)
  )
  for (i in seq_along(warn_conds)) {
    caught <- tryCatch(
      warning(warn_conds[[i]]),
      lwdid_warning = function(c) c
    )
    expect_true(inherits(caught, "lwdid_warning"),
                info = paste("Warning condition", i, "not caught by lwdid_warning"))
  }
})

# ============================================================================
# 第6组：附加信息字段测试
# ============================================================================

test_that("lwdid_missing_column_error stores column and available fields", {
  cond <- lwdid_missing_column_error("y", c("x", "z", "w"))
  expect_equal(cond$column, "y")
  expect_equal(cond$available, c("x", "z", "w"))
})

test_that("lwdid_no_treated_error stores treat_var and n_units fields", {
  cond <- lwdid_no_treated_error("d", 100)
  expect_equal(cond$treat_var, "d")
  expect_equal(cond$n_units, 100)
})

test_that("lwdid_unbalanced_panel_error stores all 4 fields", {
  cond <- lwdid_unbalanced_panel_error(5, 10, 3, 30.0)
  expect_equal(cond$min_obs, 5)
  expect_equal(cond$max_obs, 10)
  expect_equal(cond$n_incomplete_units, 3)
  expect_equal(cond$pct_unbalanced, 30.0)
})

test_that("lwdid_insufficient_pre_periods_error stores all fields including optional", {
  cond <- lwdid_insufficient_pre_periods_error(1, 2, "demean",
                                                cohort = 2005,
                                                excluded = 3)
  expect_equal(cond$available, 1)
  expect_equal(cond$required, 2)
  expect_equal(cond$rolling, "demean")
  expect_equal(cond$cohort, 2005)
  expect_equal(cond$excluded, 3)
})

test_that("lwdid_insufficient_pre_periods_error handles NULL optional fields", {
  cond <- lwdid_insufficient_pre_periods_error(1, 2, "demean")
  expect_null(cond$cohort)
  expect_null(cond$excluded)
})

test_that("lwdid_overlap_warning stores ps_range vector field", {
  cond <- lwdid_overlap_warning(5, c(0.01, 0.99), 50.0)
  expect_equal(cond$ps_range, c(0.01, 0.99))
  expect_equal(length(cond$ps_range), 2L)
})

test_that("lwdid_convergence_warning stores iteration fields", {
  cond <- lwdid_convergence_warning("logistic", 25, 100, 1e-8)
  expect_equal(cond$model, "logistic")
  expect_equal(cond$iterations, 25)
  expect_equal(cond$max_iterations, 100)
  expect_equal(cond$tolerance, 1e-8)
})

test_that("lwdid_small_sample_warning stores optional cohort/period fields", {
  cond <- lwdid_small_sample_warning(20, 8, 12, cohort = 2005, period = 2010)
  expect_equal(cond$cohort, 2005)
  expect_equal(cond$period, 2010)
  # Without optional fields
  cond2 <- lwdid_small_sample_warning(20, 8, 12)
  expect_null(cond2$cohort)
  expect_null(cond2$period)
})

test_that("lwdid_numerical_warning stores optional condition_number field", {
  cond <- lwdid_numerical_warning("issue", "func", condition_number = 1e12)
  expect_equal(cond$condition_number, 1e12)
  # Without optional field
  cond2 <- lwdid_numerical_warning("issue", "func")
  expect_null(cond2$condition_number)
})

# ============================================================================
# 第7组：conditionMessage()测试
# ============================================================================

test_that("error condition message contains key info", {
  cond <- lwdid_missing_column_error("outcome", c("x", "y"))
  msg <- conditionMessage(cond)
  expect_true(grepl("outcome", msg))
  expect_true(grepl("x, y", msg))
})

test_that("warning condition message contains key info", {
  cond <- lwdid_overlap_warning(5, c(0.01, 0.99), 50.0)
  msg <- conditionMessage(cond)
  expect_true(grepl("5 units", msg))
  expect_true(grepl("0.0100", msg))
  expect_true(grepl("50.00", msg))
})

test_that("condition message with cohort includes cohort info", {
  cond <- lwdid_small_sample_warning(20, 8, 12, cohort = 2005, period = 2010)
  msg <- conditionMessage(cond)
  expect_true(grepl("2005", msg))
  expect_true(grepl("2010", msg))
})

test_that("numerical warning with condition_number includes it in message", {
  cond <- lwdid_numerical_warning("issue", "func", condition_number = 1e12)
  msg <- conditionMessage(cond)
  expect_true(grepl("condition number", msg))
  expect_true(grepl("1.00e\\+12", msg))
})

# ============================================================================
# 第8组：消息模板精确匹配测试（验证偏差修正后的正确性）
# ============================================================================

test_that("lwdid_insufficient_data_error message starts with 'Insufficient data for estimation:'", {
  cond <- lwdid_insufficient_data_error(10, 3, 7)
  msg <- conditionMessage(cond)
  expect_true(grepl("^Insufficient data for estimation:", msg))
})

test_that("lwdid_unbalanced_panel_error pct_unbalanced is not multiplied by 100", {
  cond <- lwdid_unbalanced_panel_error(5, 10, 3, 30.0)
  msg <- conditionMessage(cond)
  expect_true(grepl("30.0%", msg, fixed = TRUE))
  expect_false(grepl("3000.0%", msg, fixed = TRUE))
})

test_that("lwdid_time_discontinuity_error message contains '(unit var:' and 'gaps at'", {
  cond <- lwdid_time_discontinuity_error(c(2005, 2008), "state", "year")
  msg <- conditionMessage(cond)
  expect_true(grepl("(unit var:", msg, fixed = TRUE))
  expect_true(grepl("gaps at", msg, fixed = TRUE))
})

test_that("lwdid_randomization_error_cond message contains 'failed after' and 'replications'", {
  cond <- lwdid_randomization_error_cond(500, "timeout")
  msg <- conditionMessage(cond)
  expect_true(grepl("failed after", msg, fixed = TRUE))
  expect_true(grepl("replications", msg, fixed = TRUE))
})

test_that("lwdid_visualization_error_cond message contains \"for '\" and \"' plot:\"", {
  cond <- lwdid_visualization_error_cond("event_study", "missing data")
  msg <- conditionMessage(cond)
  expect_true(grepl("for '", msg, fixed = TRUE))
  expect_true(grepl("' plot:", msg, fixed = TRUE))
})

test_that("lwdid_invalid_staggered_data_error message contains 'Invalid staggered data in' and 'Invalid values:'", {
  cond <- lwdid_invalid_staggered_data_error("gvar", c(-1, -2), "negative values")
  msg <- conditionMessage(cond)
  expect_true(grepl("Invalid staggered data in '", msg, fixed = TRUE))
  expect_true(grepl("Invalid values:", msg, fixed = TRUE))
})

test_that("lwdid_aggregation_error_cond message contains 'Aggregation failed for type'", {
  cond <- lwdid_aggregation_error_cond("overall", "no data")
  msg <- conditionMessage(cond)
  expect_true(grepl("Aggregation failed for type '", msg, fixed = TRUE))
})

test_that("lwdid_insufficient_cell_size_error message starts with 'All' and contains 'Output panel would be empty.'", {
  cond <- lwdid_insufficient_cell_size_error(10, 5, 20)
  msg <- conditionMessage(cond)
  expect_true(grepl("^All", msg))
  expect_true(grepl("Output panel would be empty.", msg, fixed = TRUE))
})

test_that("lwdid_small_sample_warning message contains '(N=20, treated=8, control=12)'", {
  cond <- lwdid_small_sample_warning(20, 8, 12)
  msg <- conditionMessage(cond)
  expect_true(grepl("(N=20, treated=8, control=12)", msg, fixed = TRUE))
})

test_that("lwdid_invalid_rolling_error message contains 'Must be one of:'", {
  cond <- lwdid_invalid_rolling_error("bad")
  msg <- conditionMessage(cond)
  expect_true(grepl("Must be one of:", msg, fixed = TRUE))
  expect_false(grepl("Allowed:", msg, fixed = TRUE))
})

test_that("lwdid_invalid_vce_error message contains 'Must be one of:'", {
  cond <- lwdid_invalid_vce_error("bad")
  msg <- conditionMessage(cond)
  expect_true(grepl("Must be one of:", msg, fixed = TRUE))
  expect_false(grepl("Allowed:", msg, fixed = TRUE))
})

test_that("lwdid_no_never_treated_error message contains 'with control_group='", {
  cond <- lwdid_no_never_treated_error("cohort", "never_treated")
  msg <- conditionMessage(cond)
  expect_true(grepl("with control_group=", msg, fixed = TRUE))
  expect_false(grepl("/ control_group=", msg, fixed = TRUE))
})

test_that("lwdid_numerical_warning message contains 'Numerical issue in func:'", {
  cond <- lwdid_numerical_warning("issue", "func")
  msg <- conditionMessage(cond)
  expect_true(grepl("Numerical issue in func:", msg, fixed = TRUE))
  expect_false(grepl("Numerical issue: issue", msg, fixed = TRUE))
})

test_that("lwdid_data_warning message is fixed format 'Data issue: ... Action taken: ...'", {
  cond <- lwdid_data_warning("issue", "action")
  msg <- conditionMessage(cond)
  expect_equal(msg, "Data issue: issue. Action taken: action")
})

test_that("lwdid_convergence_warning message contains 'Convergence failure in logistic:'", {
  cond <- lwdid_convergence_warning("logistic", 25, 100, 1e-8)
  msg <- conditionMessage(cond)
  expect_true(grepl("Convergence failure in logistic:", msg, fixed = TRUE))
  expect_false(grepl("Convergence issue:", msg, fixed = TRUE))
})

# ============================================================================
# 第9组：边界情况测试
# ============================================================================

test_that("lwdid_invalid_parameter_error handles empty allowed vector", {
  cond <- lwdid_invalid_parameter_error("p", "v", character(0))
  expect_s3_class(cond, "lwdid_invalid_parameter")
  msg <- conditionMessage(cond)
  expect_true(nchar(msg) > 0)
})

test_that("lwdid_insufficient_pre_periods_error handles NULL optional params", {
  cond <- lwdid_insufficient_pre_periods_error(1, 2, "demean")
  expect_s3_class(cond, "lwdid_insufficient_pre_periods")
  msg <- conditionMessage(cond)
  # Should not contain cohort or excluded info
  expect_false(grepl("Cohort", msg, fixed = TRUE))
  expect_false(grepl("excluded", msg, fixed = TRUE))
})

test_that("lwdid_invalid_staggered_data_error truncates invalid_values > 10 in message but keeps full in object", {
  long_vals <- 1:15
  cond <- lwdid_invalid_staggered_data_error("gvar", long_vals, "too many")
  msg <- conditionMessage(cond)
  # The 11th value (11) should NOT appear in the message
  # Message should contain values 1-10 but not 11-15
  expect_true(grepl("\\b10\\b", msg))
  expect_false(grepl("\\b11\\b", msg))
  # But the condition object should have the full vector
  expect_equal(length(cond$invalid_values), 15L)
  expect_equal(cond$invalid_values, 1:15)
})
