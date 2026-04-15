# ==============================================================================
# Test: export.R — LaTeX与CSV导出
# TC-10.3.1 to TC-10.3.76
# ==============================================================================

# Rebind to the frozen export-specific mock factories so helper name
# collisions from visualization/export suites do not change these tests.
.mock_common_timing_result <- .mock_export_common_timing_result
.mock_staggered_result <- .mock_export_staggered_result

# TC-10.3.1: to_latex生成合法LaTeX
test_that("TC-10.3.1: to_latex generates valid LaTeX", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x)
  expect_true(grepl("\\\\begin\\{table\\}", tex))
  expect_true(grepl("\\\\end\\{table\\}", tex))
  expect_true(grepl("\\\\begin\\{tabular\\}", tex))
})

# TC-10.3.2: 显著性星号正确
test_that("TC-10.3.2: significance stars correct for p<0.01", {
  x <- .mock_common_timing_result(pvalue = 0.001)
  tex <- to_latex(x, stars = TRUE)
  expect_true(grepl("\\*\\*\\*", tex))
})

# TC-10.3.3: to_csv summary导出单行
test_that("TC-10.3.3: to_csv summary exports single row with att", {
  x <- .mock_common_timing_result()
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "summary")
  expect_equal(nrow(df), 1L)
  expect_true("att" %in% names(df))
})

# TC-10.3.4: to_latex_comparison多模型表格
test_that("TC-10.3.4: to_latex_comparison generates multi-model table", {
  m1 <- .mock_common_timing_result(att = 1.0, pvalue = 0.01)
  m2 <- .mock_common_timing_result(att = 2.0, pvalue = 0.04)
  tex <- to_latex_comparison(m1, m2)
  expect_true(grepl("\\\\begin\\{table\\}", tex))
  expect_true(grepl("1\\.0000", tex))
  expect_true(grepl("2\\.0000", tex))
})

# TC-10.3.5: file参数写入文件
test_that("TC-10.3.5: file parameter writes to file", {
  x <- .mock_common_timing_result()
  f <- tempfile(fileext = ".tex")
  on.exit(unlink(f), add = TRUE)
  to_latex(x, file = f)
  expect_true(file.exists(f))
  content <- readLines(f)
  expect_true(any(grepl("\\\\begin\\{table\\}", content)))
})

# TC-10.3.6: to_latex_comparison每行正确以\\结尾
test_that("TC-10.3.6: comparison rows end with \\\\", {
  m1 <- .mock_common_timing_result()
  m2 <- .mock_common_timing_result(att = 2.0)
  tex <- to_latex_comparison(m1, m2)
  tex_lines <- strsplit(tex, "\n")[[1]]
  # ATT row should end with \\
  att_line <- tex_lines[grepl("^ATT", tex_lines)]
  expect_true(grepl("\\\\\\\\$", att_line))
})

# TC-10.3.7: to_latex转义下划线
test_that("TC-10.3.7: to_latex escapes underscore in estimator", {
  x <- .mock_common_timing_result(estimator = "DR_improved")
  tex <- to_latex(x)
  expect_true(grepl("DR\\\\_improved", tex))
})

# TC-10.3.8: to_latex_comparison使用命名模型
test_that("TC-10.3.8: comparison with named models", {
  m1 <- .mock_common_timing_result()
  m2 <- .mock_common_timing_result(att = 2.0)
  tex <- to_latex_comparison(m1, m2, model_names = c("Base", "Extended"))
  expect_true(grepl("Base", tex))
  expect_true(grepl("Extended", tex))
})

# TC-10.3.9: include_ci=FALSE不包含CI
test_that("TC-10.3.9: include_ci=FALSE excludes CI row", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x, include_ci = FALSE)
  expect_false(grepl("CI \\[", tex))
})

# TC-10.3.10: to_csv all导出staggered最细粒度
test_that("TC-10.3.10: to_csv all exports finest granularity", {
  x <- .mock_staggered_result()
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "all")
  # att_by_cohort_time has 9 rows (3 cohorts x 3 times)
  expect_equal(nrow(df), 9L)
})

# TC-10.3.11: stars=TRUE使用threeparttable
test_that("TC-10.3.11: stars=TRUE uses threeparttable", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x, stars = TRUE)
  expect_true(grepl("threeparttable", tex))
  expect_true(grepl("tablenotes", tex))
})

# TC-10.3.12: include_ri=TRUE输出RI信息
test_that("TC-10.3.12: include_ri=TRUE outputs RI info", {
  x <- .mock_common_timing_result(
    ri_pvalue = 0.023, ri_seed = 42L, rireps = 999L
  )
  tex <- to_latex(x, include_ri = TRUE)
  expect_true(grepl("RI p-value", tex))
  expect_true(grepl("RI seed", tex))
  expect_true(grepl("RI reps", tex))
  expect_true(grepl("42", tex))
  expect_true(grepl("999", tex))
})

# TC-10.3.13: include_ri=FALSE不输出RI信息
test_that("TC-10.3.13: include_ri=FALSE excludes RI info", {
  x <- .mock_common_timing_result(
    ri_pvalue = 0.023, ri_seed = 42L, rireps = 999L
  )
  tex <- to_latex(x, include_ri = FALSE)
  expect_false(grepl("RI p-value", tex))
})

# TC-10.3.14: include_ri=TRUE但ri_pvalue为NULL
test_that("TC-10.3.14: include_ri=TRUE with NULL ri_pvalue", {
  x <- .mock_common_timing_result(ri_pvalue = NULL)
  tex <- to_latex(x, include_ri = TRUE)
  expect_false(grepl("RI p-value", tex))
})

# TC-10.3.15: RI信息数字格式正确
test_that("TC-10.3.15: RI info number format correct", {
  x <- .mock_common_timing_result(
    ri_pvalue = 0.0234, ri_seed = 42L, rireps = 999L
  )
  tex <- to_latex(x, include_ri = TRUE, digits = 4L)
  # ri_pvalue should be formatted with 4 decimal places
  expect_true(grepl("0\\.0234", tex))
  # seed and reps should be integers
  expect_true(grepl("42 \\\\\\\\", tex))
  expect_true(grepl("999 \\\\\\\\", tex))
})

# TC-10.3.16: include_periods=TRUE输出时期效应表格
test_that("TC-10.3.16: include_periods=TRUE outputs period table", {
  abp <- data.frame(
    period = c(-2L, -1L, 0L),
    att = c(0.1, 0.2, 1.5),
    se = c(0.05, 0.06, 0.3),
    t_stat = c(2.0, 3.33, 5.0),
    pvalue = c(0.05, 0.001, 0.0001),
    ci_lower = c(0.0, 0.08, 0.9),
    ci_upper = c(0.2, 0.32, 2.1),
    stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(att_by_period = abp)
  tex <- to_latex(x, include_periods = TRUE)
  expect_true(grepl("Period & ATT & SE & t-stat & p-value & CI", tex))
  expect_true(grepl("-2", tex))
})

# TC-10.3.17: include_periods=FALSE不输出时期效应表格
test_that("TC-10.3.17: include_periods=FALSE excludes period table", {
  abp <- data.frame(
    period = c(-2L, -1L, 0L),
    att = c(0.1, 0.2, 1.5), se = c(0.05, 0.06, 0.3),
    pvalue = c(0.05, 0.001, 0.0001),
    ci_lower = c(0.0, 0.08, 0.9), ci_upper = c(0.2, 0.32, 2.1),
    stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(att_by_period = abp)
  tex <- to_latex(x, include_periods = FALSE)
  expect_false(grepl("Period & ATT & SE", tex))
})

# TC-10.3.18: 时期效应表格显著性星号
test_that("TC-10.3.18: period table significance stars", {
  abp <- data.frame(
    period = c(-1L, 0L),
    att = c(0.2, 1.5), se = c(0.06, 0.3),
    pvalue = c(0.001, 0.04),
    ci_lower = c(0.08, 0.9), ci_upper = c(0.32, 2.1),
    stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(att_by_period = abp)
  tex <- to_latex(x, include_periods = TRUE, stars = TRUE)
  # p=0.001 < 0.01 → ***
  expect_true(grepl("\\*\\*\\*", tex))
  # p=0.04 < 0.05 → **
  expect_true(grepl("\\*\\*[^*]", tex))
})

# TC-10.3.19: 时期效应表格处理conf.low/conf.high列名变体
test_that("TC-10.3.19: period table handles conf.low/conf.high", {
  abp <- data.frame(
    event_time = c(-1L, 0L),
    att = c(0.2, 1.5), se = c(0.06, 0.3),
    pvalue = c(0.05, 0.001),
    conf.low = c(0.08, 0.9), conf.high = c(0.32, 2.1),
    stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(att_by_period = abp)
  tex <- to_latex(x, include_periods = TRUE)
  # Should contain CI values, not "--"
  expect_true(grepl("\\[0\\.0800", tex))
})

# TC-10.3.20: include_periods=TRUE但att_by_period为NULL
test_that("TC-10.3.20: include_periods=TRUE with NULL att_by_period", {
  x <- .mock_common_timing_result(att_by_period = NULL)
  tex <- to_latex(x, include_periods = TRUE)
  # Should not crash, just no period table
  expect_false(grepl("Period & ATT & SE", tex))
})

# TC-10.3.21: include_staggered=TRUE输出队列效应表格
test_that("TC-10.3.21: include_staggered=TRUE outputs cohort table", {
  x <- .mock_staggered_result()
  tex <- to_latex(x, include_staggered = TRUE)
  expect_true(grepl("Cohort & ATT & SE & N & Weight", tex))
  expect_true(grepl("Overall \\(weighted\\)", tex))
})

# TC-10.3.22: include_staggered=FALSE不输出队列效应表格
test_that("TC-10.3.22: include_staggered=FALSE excludes cohort table", {
  x <- .mock_staggered_result()
  tex <- to_latex(x, include_staggered = FALSE)
  expect_false(grepl("Cohort & ATT & SE & N & Weight", tex))
})

# TC-10.3.23: include_staggered=TRUE但非staggered结果
test_that("TC-10.3.23: include_staggered=TRUE with non-staggered", {
  x <- .mock_common_timing_result(is_staggered = FALSE)
  tex <- to_latex(x, include_staggered = TRUE)
  expect_false(grepl("Cohort & ATT & SE & N & Weight", tex))
})

# TC-10.3.24: 队列效应表格权重从cohort_weights回退获取
test_that("TC-10.3.24: cohort weights fallback to cohort_weights", {
  abc <- data.frame(
    cohort = c(2004L, 2006L),
    att = c(1.2, 1.8), se = c(0.25, 0.30),
    n = c(100L, 150L),
    stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(
    is_staggered = TRUE, att_by_cohort = abc,
    cohort_weights = c("2004" = 0.40, "2006" = 0.60)
  )
  tex <- to_latex(x, include_staggered = TRUE)
  expect_true(grepl("0\\.4000", tex))
  expect_true(grepl("0\\.6000", tex))
})

# TC-10.3.25: 队列效应表格Overall行ATT与主表格ATT一致
test_that("TC-10.3.25: cohort Overall ATT matches main ATT", {
  x <- .mock_staggered_result()
  tex <- to_latex(x, include_staggered = TRUE, digits = 4L)
  # Overall row should contain the main ATT (1.5000)
  overall_match <- regmatches(tex, regexpr("Overall \\(weighted\\)[^\n]+", tex))
  expect_true(grepl("1\\.5000", overall_match))
})

# TC-10.3.26: 同时启用include_ri+include_periods+include_staggered
test_that("TC-10.3.26: all optional sections enabled", {
  abp <- data.frame(
    period = c(-1L, 0L), att = c(0.2, 1.5), se = c(0.06, 0.3),
    pvalue = c(0.05, 0.001),
    ci_lower = c(0.08, 0.9), ci_upper = c(0.32, 2.1),
    stringsAsFactors = FALSE
  )
  x <- .mock_staggered_result()
  x$ri_pvalue <- 0.023
  x$ri_seed <- 42L
  x$rireps <- 999L
  x$att_by_period <- abp
  tex <- to_latex(x, include_ri = TRUE, include_periods = TRUE,
                  include_staggered = TRUE)
  expect_true(grepl("RI p-value", tex))
  expect_true(grepl("Period & ATT", tex))
  expect_true(grepl("Cohort & ATT", tex))
})

# TC-10.3.27: 默认参数向后兼容
test_that("TC-10.3.27: default params backward compatible", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x)
  # include_ri=TRUE by default, but no ri_pvalue → no RI section
  expect_false(grepl("RI p-value", tex))
  # include_periods=FALSE by default
  expect_false(grepl("Period & ATT", tex))
  # include_staggered=FALSE by default
  expect_false(grepl("Cohort & ATT", tex))
})

# TC-10.3.28: .escape_latex转义所有特殊字符
test_that("TC-10.3.28: .escape_latex escapes all special chars", {
  esc <- lwdid:::.escape_latex
  expect_equal(esc("a_b"), "a\\_b")
  expect_equal(esc("a%b"), "a\\%b")
  expect_equal(esc("a&b"), "a\\&b")
  expect_equal(esc("a#b"), "a\\#b")
  expect_equal(esc("a$b"), "a\\$b")
  expect_equal(esc("a{b}"), "a\\{b\\}")
  expect_true(grepl("textasciitilde", esc("a~b")))
  expect_true(grepl("textasciicircum", esc("a^b")))
  expect_true(grepl("textbackslash", esc("a\\b")))
})

# TC-10.3.29: 时期效应表格p<0.001显示为$<$0.001
test_that("TC-10.3.29: period table p<0.001 format", {
  abp <- data.frame(
    period = 0L, att = 1.5, se = 0.3, pvalue = 0.0001,
    ci_lower = 0.9, ci_upper = 2.1, stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(att_by_period = abp)
  tex <- to_latex(x, include_periods = TRUE)
  expect_true(grepl("\\$<\\$0\\.001", tex))
})

# TC-10.3.30: to_csv by_cohort在非staggered时报错
test_that("TC-10.3.30: to_csv by_cohort errors for non-staggered", {
  x <- .mock_common_timing_result()
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  expect_error(to_csv(x, file = f, what = "by_cohort"),
               "No cohort-specific results")
})

# TC-10.3.31: to_csv summary包含RI信息
test_that("TC-10.3.31: to_csv summary includes RI columns", {
  x <- .mock_common_timing_result(
    ri_pvalue = 0.023, ri_seed = 42L, rireps = 999L
  )
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "summary")
  expect_true("ri_pvalue" %in% names(df))
  expect_true("ri_seed" %in% names(df))
  expect_true("rireps" %in% names(df))
  expect_equal(df$ri_pvalue, 0.023)
})

# TC-10.3.32: to_csv summary无RI时不包含RI列
test_that("TC-10.3.32: to_csv summary without RI has no RI columns", {
  x <- .mock_common_timing_result(ri_pvalue = NULL)
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "summary")
  expect_false("ri_pvalue" %in% names(df))
})

# TC-10.3.33: to_latex_comparison转义model_names特殊字符
test_that("TC-10.3.33: comparison escapes model_names", {
  m1 <- .mock_common_timing_result()
  m2 <- .mock_common_timing_result(att = 2.0)
  tex <- to_latex_comparison(m1, m2,
                              model_names = c("Model_A", "Model_B"))
  expect_true(grepl("Model\\\\_A", tex))
  expect_true(grepl("Model\\\\_B", tex))
})

# TC-10.3.34: to_latex_comparison包含VCE行
test_that("TC-10.3.34: comparison includes VCE row", {
  m1 <- .mock_common_timing_result(vce_type = "HC1_robust")
  m2 <- .mock_common_timing_result(vce_type = "cluster")
  tex <- to_latex_comparison(m1, m2)
  expect_true(grepl("VCE", tex))
  expect_true(grepl("HC1\\\\_robust", tex))
  expect_true(grepl("cluster", tex))
})

# TC-10.3.35: to_latex digits参数验证
test_that("TC-10.3.35: to_latex digits validation", {
  x <- .mock_common_timing_result()
  expect_error(to_latex(x, digits = -1))
  expect_error(to_latex(x, digits = "abc"))
})

# TC-10.3.36: to_latex逻辑参数类型验证
test_that("TC-10.3.36: to_latex logical param validation", {
  x <- .mock_common_timing_result()
  expect_error(to_latex(x, include_ci = "yes"))
  expect_error(to_latex(x, booktabs = 1))
  expect_error(to_latex(x, stars = "true"))
})

# TC-10.3.37: to_latex_comparison stars=TRUE输出星号
test_that("TC-10.3.37: comparison stars=TRUE outputs stars", {
  m1 <- .mock_common_timing_result(pvalue = 0.001)
  m2 <- .mock_common_timing_result(pvalue = 0.04)
  tex <- to_latex_comparison(m1, m2, stars = TRUE)
  expect_true(grepl("threeparttable", tex))
  expect_true(grepl("tablenotes", tex))
  expect_true(grepl("\\*\\*\\*", tex))
})

# TC-10.3.38: to_latex_comparison stars=FALSE不输出星号
test_that("TC-10.3.38: comparison stars=FALSE no stars", {
  m1 <- .mock_common_timing_result(pvalue = 0.001)
  tex <- to_latex_comparison(m1, stars = FALSE)
  expect_false(grepl("threeparttable", tex))
  expect_false(grepl("tablenotes", tex))
  expect_false(grepl("\\*\\*\\*", tex))
})

# TC-10.3.39: to_latex_comparison include_ci=TRUE输出CI行
test_that("TC-10.3.39: comparison include_ci=TRUE outputs CI", {
  m1 <- .mock_common_timing_result()
  tex <- to_latex_comparison(m1, include_ci = TRUE)
  expect_true(grepl("CI \\[", tex))
})

# TC-10.3.40: to_latex_comparison include_ci=FALSE不输出CI行
test_that("TC-10.3.40: comparison include_ci=FALSE excludes CI", {
  m1 <- .mock_common_timing_result()
  tex <- to_latex_comparison(m1, include_ci = FALSE)
  expect_false(grepl("CI \\[", tex))
})

# TC-10.3.41: booktabs=FALSE使用hline
test_that("TC-10.3.41: booktabs=FALSE uses hline", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x, booktabs = FALSE)
  expect_true(grepl("\\\\hline", tex))
  expect_false(grepl("\\\\toprule", tex))
  expect_false(grepl("\\\\midrule", tex))
  expect_false(grepl("\\\\bottomrule", tex))
})

# TC-10.3.42: stars=FALSE不输出星号和tablenotes
test_that("TC-10.3.42: stars=FALSE no tablenotes", {
  x <- .mock_common_timing_result(pvalue = 0.001)
  tex <- to_latex(x, stars = FALSE)
  expect_false(grepl("\\*\\*\\*", tex))
  expect_false(grepl("tablenotes", tex))
})

# TC-10.3.43: se_in_parentheses=FALSE不使用括号
test_that("TC-10.3.43: se_in_parentheses=FALSE no parens", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x, se_in_parentheses = FALSE)
  expect_true(grepl("^SE &", strsplit(tex, "\n")[[1]][
    grepl("^SE &", strsplit(tex, "\n")[[1]])
  ]))
})

# TC-10.3.44: caption和label正确输出
test_that("TC-10.3.44: caption and label output", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x, caption = "My Table", label = "tab:main")
  expect_true(grepl("\\\\caption\\{My Table\\}", tex))
  expect_true(grepl("\\\\label\\{tab:main\\}", tex))
})

# TC-10.3.45: to_csv by_period正常导出
test_that("TC-10.3.45: to_csv by_period exports correctly", {
  abp <- data.frame(
    period = c(-2L, -1L, 0L),
    att = c(0.1, 0.2, 1.5), se = c(0.05, 0.06, 0.3),
    pvalue = c(0.05, 0.001, 0.0001),
    stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(att_by_period = abp)
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "by_period")
  expect_equal(nrow(df), 3L)
  expect_true("att" %in% names(df))
  expect_true("period" %in% names(df))
})

# TC-10.3.46: to_csv what参数非法值报错
test_that("TC-10.3.46: to_csv invalid what errors", {
  x <- .mock_common_timing_result()
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  expect_error(to_csv(x, file = f, what = "invalid"))
})

# TC-10.3.47: to_csv summary包含alpha字段
test_that("TC-10.3.47: to_csv summary includes alpha", {
  x <- .mock_common_timing_result(alpha = 0.05)
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "summary")
  expect_true("alpha" %in% names(df))
  expect_equal(df$alpha, 0.05)
})

# TC-10.3.48: to_latex_comparison单模型正确生成
test_that("TC-10.3.48: comparison single model generates {lc}", {
  m1 <- .mock_common_timing_result()
  tex <- to_latex_comparison(m1)
  expect_true(grepl("\\{lc\\}", tex))
})

# TC-10.3.49: to_csv by_period att_by_period为NULL时报错
test_that("TC-10.3.49: to_csv by_period NULL att_by_period errors", {
  x <- .mock_common_timing_result(att_by_period = NULL)
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  expect_error(to_csv(x, file = f, what = "by_period"))
})

# TC-10.3.50: to_csv by_period空data.frame时报错
test_that("TC-10.3.50: to_csv by_period empty df errors", {
  x <- .mock_common_timing_result(
    att_by_period = data.frame(period = integer(0), att = numeric(0),
                                se = numeric(0))
  )
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  expect_error(to_csv(x, file = f, what = "by_period"))
})

# TC-10.3.51: to_csv all所有数据源为NULL时回退
test_that("TC-10.3.51: to_csv all fallback to minimal df", {
  x <- .mock_common_timing_result(
    att_by_cohort_time = NULL, att_by_period = NULL
  )
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "all")
  expect_equal(nrow(df), 1L)
  expect_true("att" %in% names(df))
  expect_true("se" %in% names(df))
})

# TC-10.3.52: rolling逻辑值正确处理
test_that("TC-10.3.52: rolling logical handled in latex and csv", {
  x <- .mock_common_timing_result(rolling = TRUE)
  tex <- to_latex(x)
  expect_true(grepl("TRUE", tex))
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "summary")
  expect_equal(df$rolling, "TRUE")
})

# TC-10.3.53: to_csv by_cohort空data.frame时报错
test_that("TC-10.3.53: to_csv by_cohort empty df errors", {
  x <- .mock_common_timing_result(
    att_by_cohort = data.frame(cohort = integer(0), att = numeric(0),
                                se = numeric(0))
  )
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  expect_error(to_csv(x, file = f, what = "by_cohort"),
               "No cohort-specific results")
})

# TC-10.3.54: to_latex_comparison参数验证
test_that("TC-10.3.54: comparison param validation", {
  m1 <- .mock_common_timing_result()
  expect_error(to_latex_comparison(m1, digits = -1))
  expect_error(to_latex_comparison(m1, include_ci = "yes"))
})

# TC-10.3.55: to_latex_comparison 0个模型报错
test_that("TC-10.3.55: comparison zero models errors", {
  expect_error(to_latex_comparison())
})

# TC-10.3.56: to_latex_comparison model_names长度不匹配报错
test_that("TC-10.3.56: comparison model_names length mismatch", {
  m1 <- .mock_common_timing_result()
  m2 <- .mock_common_timing_result()
  expect_error(to_latex_comparison(m1, m2, model_names = c("A")))
})

# TC-10.3.57: to_latex_comparison传入非lwdid_result报错
test_that("TC-10.3.57: comparison non-lwdid_result errors", {
  expect_error(to_latex_comparison(list(a = 1)))
})

# TC-10.3.58: include_diagnostics=TRUE输出诊断信息
test_that("TC-10.3.58: include_diagnostics=TRUE outputs diagnostics", {
  diag <- list(
    trends = list(recommended_method = "demean"),
    selection = list(selection_risk = "Low")
  )
  x <- .mock_common_timing_result(diagnostics = diag)
  tex <- to_latex(x, include_diagnostics = TRUE)
  expect_true(grepl("Recommendation", tex))
  expect_true(grepl("demean", tex))
  expect_true(grepl("Selection Risk", tex))
  expect_true(grepl("Low", tex))
})

# TC-10.3.59: include_diagnostics=TRUE但diagnostics为NULL
test_that("TC-10.3.59: include_diagnostics=TRUE with NULL diagnostics", {
  x <- .mock_common_timing_result(diagnostics = NULL)
  tex <- to_latex(x, include_diagnostics = TRUE)
  expect_false(grepl("Recommendation", tex))
})

# TC-10.3.60: se_in_parentheses=TRUE时SE括号格式
test_that("TC-10.3.60: se_in_parentheses=TRUE shows parens", {
  x <- .mock_common_timing_result(se_att = 0.3)
  tex <- to_latex(x, se_in_parentheses = TRUE, digits = 4L)
  expect_true(grepl("\\(0\\.3000\\)", tex))
})

# TC-10.3.61: booktabs=TRUE使用toprule/midrule/bottomrule
test_that("TC-10.3.61: booktabs=TRUE uses booktabs rules", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x, booktabs = TRUE)
  expect_true(grepl("\\\\toprule", tex))
  expect_true(grepl("\\\\midrule", tex))
  expect_true(grepl("\\\\bottomrule", tex))
})

# TC-10.3.62: comparison文本字段转义
test_that("TC-10.3.62: comparison escapes text fields", {
  m1 <- .mock_common_timing_result(
    estimator = "DR_v2", rolling = TRUE, vce_type = "HC1_robust"
  )
  tex <- to_latex_comparison(m1)
  expect_true(grepl("DR\\\\_v2", tex))
  expect_true(grepl("HC1\\\\_robust", tex))
  expect_true(grepl("TRUE", tex))
})

# TC-10.3.63: float参数验证
test_that("TC-10.3.63: float parameter validation", {
  x <- .mock_common_timing_result()
  # Empty string is valid character
  tex <- to_latex(x, float = "")
  expect_true(grepl("\\\\begin\\{table\\}\\[\\]", tex))
})

# TC-10.3.64: pvalue为NA时星号为空且p-value显示--
test_that("TC-10.3.64: NA pvalue shows -- and no stars", {
  x <- .mock_common_timing_result(pvalue = NA_real_)
  tex <- to_latex(x, stars = TRUE)
  # No stars on ATT line
  att_line <- strsplit(tex, "\n")[[1]][grepl("^ATT", strsplit(tex, "\n")[[1]])]
  expect_false(grepl("\\*", att_line))
  # p-value row shows --
  expect_true(grepl("p-value & --", tex))
})

# TC-10.3.65: comparison中某模型pvalue为NA
test_that("TC-10.3.65: comparison NA pvalue model", {
  m1 <- .mock_common_timing_result(pvalue = 0.001)
  m2 <- .mock_common_timing_result(pvalue = NA_real_)
  tex <- to_latex_comparison(m1, m2, stars = TRUE)
  # m1 should have ***, m2 should have no stars
  expect_true(grepl("\\*\\*\\*", tex))
  # p-value row should contain --
  pval_line <- strsplit(tex, "\n")[[1]][grepl("^p-value", strsplit(tex, "\n")[[1]])]
  expect_true(grepl("--", pval_line))
})

# TC-10.3.66: 时期效应表格pvalue为NA
test_that("TC-10.3.66: period table NA pvalue shows --", {
  abp <- data.frame(
    period = c(-1L, 0L),
    att = c(0.2, 1.5), se = c(0.06, 0.3),
    pvalue = c(NA_real_, 0.0001),
    ci_lower = c(0.08, 0.9), ci_upper = c(0.32, 2.1),
    stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(att_by_period = abp)
  tex <- to_latex(x, include_periods = TRUE, stars = TRUE)
  tex_lines <- strsplit(tex, "\n")[[1]]
  # Find the period -1 line (NA pvalue) — should have "--" for p-value
  period_neg1 <- tex_lines[grepl("^-1 &", tex_lines)]
  expect_true(grepl("--", period_neg1))
  # Period 0 line should have $<$0.001 (p=0.0001 < 0.001)
  period_0 <- tex_lines[grepl("^0 &", tex_lines)]
  expect_true(grepl("\\$<\\$0\\.001", period_0))
})

# TC-10.3.67: .escape_latex接受数值输入
test_that("TC-10.3.67: .escape_latex accepts numeric input", {
  esc <- lwdid:::.escape_latex
  expect_equal(esc(123), "123")
  expect_equal(esc(0.456), "0.456")
})

# TC-10.3.68: .escape_latex接受逻辑值输入
test_that("TC-10.3.68: .escape_latex accepts logical input", {
  esc <- lwdid:::.escape_latex
  expect_equal(esc(TRUE), "TRUE")
  expect_equal(esc(FALSE), "FALSE")
})

# TC-10.3.69: to_csv summary rolling=TRUE为字符串
test_that("TC-10.3.69: to_csv rolling=TRUE is string not integer", {
  x <- .mock_common_timing_result(rolling = TRUE)
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "summary")
  expect_true(is.character(df$rolling))
  expect_equal(df$rolling, "TRUE")
  # Verify in actual CSV file
  csv_content <- readLines(f)
  # rolling column should contain "TRUE" not "1"
  expect_true(any(grepl("TRUE", csv_content)))
})

# TC-10.3.70: 时期效应子表格booktabs=FALSE使用hline
test_that("TC-10.3.70: period subtable booktabs=FALSE uses hline", {
  abp <- data.frame(
    period = 0L, att = 1.5, se = 0.3, pvalue = 0.001,
    ci_lower = 0.9, ci_upper = 2.1, stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(att_by_period = abp)
  tex <- to_latex(x, include_periods = TRUE, booktabs = FALSE)
  # Should use \hline throughout, no \toprule
  expect_false(grepl("\\\\toprule", tex))
  expect_true(grepl("\\\\hline", tex))
})

# TC-10.3.71: 队列效应子表格booktabs=FALSE使用hline
test_that("TC-10.3.71: cohort subtable booktabs=FALSE uses hline", {
  x <- .mock_staggered_result()
  tex <- to_latex(x, include_staggered = TRUE, booktabs = FALSE)
  expect_false(grepl("\\\\toprule", tex))
  expect_true(grepl("\\\\hline", tex))
})

# TC-10.3.72: to_csv file=NULL报错
test_that("TC-10.3.72: to_csv file=NULL errors", {
  x <- .mock_common_timing_result()
  expect_error(to_csv(x, file = NULL))
})

# TC-10.3.73: to_csv all回退到att_by_period
test_that("TC-10.3.73: to_csv all fallback to att_by_period", {
  abp <- data.frame(
    period = c(-1L, 0L), att = c(0.2, 1.5), se = c(0.06, 0.3),
    stringsAsFactors = FALSE
  )
  x <- .mock_common_timing_result(
    att_by_cohort_time = NULL, att_by_period = abp
  )
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  df <- to_csv(x, file = f, what = "all")
  expect_equal(nrow(df), 2L)
  expect_true("period" %in% names(df))
})

# TC-10.3.74: comparison SE始终括号格式
test_that("TC-10.3.74: comparison SE always in parentheses", {
  m1 <- .mock_common_timing_result(se_att = 0.3)
  tex <- to_latex_comparison(m1)
  expect_true(grepl("\\(0\\.3000\\)", tex))
})

# TC-10.3.75: to_latex stars=FALSE仍包含threeparttable
test_that("TC-10.3.75: to_latex stars=FALSE still has threeparttable", {
  x <- .mock_common_timing_result()
  tex <- to_latex(x, stars = FALSE)
  expect_true(grepl("threeparttable", tex))
  expect_false(grepl("tablenotes", tex))
})

# TC-10.3.76: to_latex_comparison stars=FALSE不包含threeparttable
test_that("TC-10.3.76: comparison stars=FALSE no threeparttable", {
  m1 <- .mock_common_timing_result()
  tex <- to_latex_comparison(m1, stars = FALSE)
  expect_false(grepl("threeparttable", tex))
  expect_false(grepl("tablenotes", tex))
})
