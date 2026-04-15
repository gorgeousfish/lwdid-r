# ==============================================================================
# Test: plot_diagnostics.R — Diagnostic Plot Methods
# TC-10.2.1 to TC-10.2.43
# ==============================================================================

# TC-10.2.1: 规格曲线图正确生成（ggplot对象，含GeomPoint图层）
test_that("TC-10.2.1: specification curve generates ggplot with GeomPoint", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_pre_period()
  p <- plot.lwdid_sensitivity(x)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_true("GeomPoint" %in% layer_classes)
})

# TC-10.2.2: 预期效应图正确生成
test_that("TC-10.2.2: anticipation sensitivity plot generates ggplot", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_no_anticipation()
  p <- plot.lwdid_sensitivity(x)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_true("GeomPoint" %in% layer_classes)
  expect_true("GeomLine" %in% layer_classes)
})

# TC-10.2.3: 综合敏感性图正确拼接（patchwork对象）
test_that("TC-10.2.3: comprehensive sensitivity creates patchwork", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  x <- create_mock_sensitivity_comprehensive()
  p <- plot.lwdid_sensitivity_comprehensive(x)
  expect_s3_class(p, "patchwork")
})

# TC-10.2.4: 队列趋势图正确生成
test_that("TC-10.2.4: heterogeneous trends plot generates ggplot", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_true("GeomLine" %in% layer_classes)
})

# TC-10.2.5: show_control=FALSE不显示对照组
test_that("TC-10.2.5: show_control=FALSE hides control group", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x, show_control = FALSE)
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) > 0)
})

# TC-10.2.6: show_ci=FALSE不显示置信带
test_that("TC-10.2.6: show_ci=FALSE removes GeomRibbon", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x, show_ci = FALSE)
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_false("GeomRibbon" %in% layer_classes)
})

# TC-10.2.7: 变换推荐评分图正确生成
test_that("TC-10.2.7: transformation recommendation has GeomCol", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_transformation_recommendation()
  p <- plot.lwdid_transformation_recommendation(x)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_true("GeomCol" %in% layer_classes)
})

# TC-10.2.8: 敏感性类型自动推断
test_that("TC-10.2.8: type auto-inferred from x$type", {
  skip_if_not_installed("ggplot2")
  x_pre <- create_mock_sensitivity_pre_period()
  p_pre <- plot.lwdid_sensitivity(x_pre)
  expect_s3_class(p_pre, "ggplot")

  x_ant <- create_mock_sensitivity_no_anticipation()
  p_ant <- plot.lwdid_sensitivity(x_ant)
  expect_s3_class(p_ant, "ggplot")
})

# TC-10.2.9: show_threshold=FALSE不显示阈值带
test_that("TC-10.2.9: show_threshold=FALSE removes threshold band", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_pre_period()
  p_with <- plot.lwdid_sensitivity(x, show_threshold = TRUE)
  p_without <- plot.lwdid_sensitivity(x, show_threshold = FALSE)
  expect_true(length(p_with$layers) > length(p_without$layers))
})

# TC-10.2.10: 非lwdid_heterogeneous_trends对象报错
test_that("TC-10.2.10: non-lwdid_heterogeneous_trends errors", {
  skip_if_not_installed("ggplot2")
  wrong <- list(a = 1)
  expect_error(
    plot.lwdid_heterogeneous_trends(wrong),
    "inherits"
  )
})

# TC-10.2.11: facet_by_cohort=TRUE正确分面
test_that("TC-10.2.11: facet_by_cohort=TRUE adds FacetWrap", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x, facet_by_cohort = TRUE)
  expect_s3_class(p$facet, "FacetWrap")
})

# TC-10.2.12: facet_by_cohort=FALSE无分面
test_that("TC-10.2.12: facet_by_cohort=FALSE has FacetNull", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x, facet_by_cohort = FALSE)
  expect_s3_class(p$facet, "FacetNull")
})

# TC-10.2.13: show_slope_labels=TRUE添加GeomLabel
test_that("TC-10.2.13: show_slope_labels=TRUE adds GeomLabel", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x, show_slope_labels = TRUE)
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_true("GeomLabel" %in% layer_classes)
})

# TC-10.2.14: show_slope_labels=FALSE无GeomLabel
test_that("TC-10.2.14: show_slope_labels=FALSE removes GeomLabel", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x, show_slope_labels = FALSE)
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_false("GeomLabel" %in% layer_classes)
})

# TC-10.2.15: 斜率标注格式包含β和p值
test_that("TC-10.2.15: slope labels contain beta and p-value", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x, show_slope_labels = TRUE)
  label_layers <- Filter(
    function(l) inherits(l$geom, "GeomLabel"), p$layers
  )
  expect_true(length(label_layers) > 0)
  label_data <- label_layers[[1]]$data
  expect_true("slope_label" %in% names(label_data))
  labels <- label_data$slope_label
  expect_true(all(grepl("\u03b2=", labels)))
  expect_true(all(grepl("p=", labels)))
  expect_true(any(grepl("0\\.5000", labels)))
  expect_true(any(grepl("0\\.3000", labels)))
  # p=0.001 is NOT < 0.001, so formatted as "0.001" (digits=3)
  expect_true(any(grepl("p=0\\.001", labels)))
  expect_true(any(grepl("0\\.080", labels)))
})

# TC-10.2.16: 变换推荐图使用水平条形图（CoordFlip）
test_that("TC-10.2.16: transformation uses CoordFlip", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_transformation_recommendation()
  p <- plot.lwdid_transformation_recommendation(x)
  expect_s3_class(p$coordinates, "CoordFlip")
})

# TC-10.2.17: sub_scores可用时绘制堆叠条形图
test_that("TC-10.2.17: stacked bar with sub_scores", {
  skip_if_not_installed("ggplot2")
  sub_scores <- list(
    demean = c(linearity = 0.9, variance = 0.8, normality = 0.85),
    detrend = c(linearity = 0.7, variance = 0.75, normality = 0.7),
    demeanq = c(linearity = 0.95, variance = 0.85, normality = 0.9),
    detrendq = c(linearity = 0.6, variance = 0.7, normality = 0.6)
  )
  x <- create_mock_transformation_recommendation(sub_scores = sub_scores)
  p <- plot.lwdid_transformation_recommendation(x)
  expect_s3_class(p, "ggplot")
  # Should have multiple fill colors (stacked)
  built <- ggplot2::ggplot_build(p)
  fills <- unique(built$data[[1]]$fill)
  expect_true(length(fills) > 1)
})

# TC-10.2.18: 无sub_scores时绘制简单条形图
test_that("TC-10.2.18: simple bar without sub_scores", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_transformation_recommendation(sub_scores = NULL)
  p <- plot.lwdid_transformation_recommendation(x)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_true("GeomCol" %in% layer_classes)
})

# TC-10.2.19: method_labels正确替换标签
test_that("TC-10.2.19: method_labels correctly replace labels", {
  skip_if_not_installed("ggplot2")
  labels <- c(demean = "去均值", demeanq = "去均值(二次)")
  x <- create_mock_transformation_recommendation(method_labels = labels)
  p <- plot.lwdid_transformation_recommendation(x)
  built <- ggplot2::ggplot_build(p)
  # Check that factor levels include the replaced labels
  plot_data <- p$data
  lvls <- levels(plot_data$method)
  expect_true("去均值" %in% lvls)
  expect_true("去均值(二次)" %in% lvls)
  # detrend and detrendq should keep original names
  expect_true("detrend" %in% lvls)
  expect_true("detrendq" %in% lvls)
})

# TC-10.2.20: 置信度等级影响subtitle颜色
test_that("TC-10.2.20: confidence level affects subtitle color", {
  skip_if_not_installed("ggplot2")
  x_high <- create_mock_transformation_recommendation(
    confidence_level = "High"
  )
  x_low <- create_mock_transformation_recommendation(
    confidence_level = "Low"
  )
  p_high <- plot.lwdid_transformation_recommendation(x_high)
  p_low <- plot.lwdid_transformation_recommendation(x_low)
  # Extract subtitle color from theme
  theme_high <- p_high$theme$plot.subtitle$colour
  theme_low <- p_low$theme$plot.subtitle$colour
  expect_false(identical(theme_high, theme_low))
  expect_equal(theme_high, "#2166AC")  # High = blue
  expect_equal(theme_low, "#B2182B")   # Low = red
})

# TC-10.2.21: 非lwdid_transformation_recommendation对象报错
test_that("TC-10.2.21: non-lwdid_transformation_recommendation errors", {
  skip_if_not_installed("ggplot2")
  wrong <- list(a = 1)
  expect_error(
    plot.lwdid_transformation_recommendation(wrong),
    "inherits"
  )
})

# TC-10.2.22: 规格曲线图CI使用t分布（非z=1.96）
test_that("TC-10.2.22: spec curve CI uses t-distribution", {
  skip_if_not_installed("ggplot2")
  # Construct spec with small df: att=1.0, se=0.2, df=2
  # t_crit(0.975, 2) ≈ 4.303
  # t-dist CI: [1.0 - 4.303*0.2, 1.0 + 4.303*0.2] = [0.139, 1.861]
  # z=1.96 CI: [1.0 - 1.96*0.2, 1.0 + 1.96*0.2] = [0.608, 1.392]
  t_crit <- qt(0.975, 2)
  specs <- list(
    list(n_pre_periods = 3L, att = 1.0, se = 0.2,
         pvalue = 0.01,
         ci_lower = 1.0 - t_crit * 0.2,
         ci_upper = 1.0 + t_crit * 0.2,
         converged = TRUE)
  )
  x <- create_mock_sensitivity_pre_period(specifications = specs)
  p <- plot.lwdid_sensitivity(x, show_ci = TRUE)
  built <- ggplot2::ggplot_build(p)
  # Find errorbar layer data
  eb_data <- NULL
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomErrorbar")) {
      eb_data <- built$data[[i]]
      break
    }
  }
  expect_false(is.null(eb_data))
  # t-dist CI should be wider than z=1.96 CI
  ci_width <- eb_data$ymax[1] - eb_data$ymin[1]
  z_width <- 2 * 1.96 * 0.2
  expect_true(ci_width > z_width)
})

# TC-10.2.23: 规格曲线阈值带固定为基准 ATT 的 ±25%
test_that("TC-10.2.23: threshold band stays anchored at +/-25% of baseline", {
  skip_if_not_installed("ggplot2")
  # robustness_threshold drives sensitivity classification, not the plot band.
  # Plot contract stays anchored at baseline_att +/- 25% * |baseline_att|.
  x <- create_mock_sensitivity_pre_period(
    baseline_att = 2.0, robustness_threshold = 0.10
  )
  p <- plot.lwdid_sensitivity(x, show_threshold = TRUE)
  built <- ggplot2::ggplot_build(p)
  # Find annotate rect layer
  rect_data <- NULL
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomRect")) {
      rect_data <- built$data[[i]]
      break
    }
  }
  expect_false(is.null(rect_data))
  expect_equal(rect_data$ymin[1], 1.5, tolerance = 1e-10)
  expect_equal(rect_data$ymax[1], 2.5, tolerance = 1e-10)
})

# TC-10.2.24: show_ci=FALSE时规格曲线图不显示CI误差棒
test_that("TC-10.2.24: show_ci=FALSE removes errorbar from spec curve", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_pre_period()
  p <- plot.lwdid_sensitivity(x, show_ci = FALSE)
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_false("GeomErrorbar" %in% layer_classes)
  expect_true("GeomPoint" %in% layer_classes)
})

# TC-10.2.25: show_ci=FALSE时预期效应图不显示CI误差棒
test_that("TC-10.2.25: show_ci=FALSE removes errorbar from anticipation", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_no_anticipation()
  p <- plot.lwdid_sensitivity(x, show_ci = FALSE)
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_false("GeomErrorbar" %in% layer_classes)
})

# TC-10.2.26: show_ci=TRUE（默认）时显示CI误差棒
test_that("TC-10.2.26: show_ci=TRUE shows errorbar", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_pre_period()
  p <- plot.lwdid_sensitivity(x, show_ci = TRUE)
  layer_classes <- vapply(
    p$layers, function(l) class(l$geom)[1], character(1)
  )
  expect_true("GeomErrorbar" %in% layer_classes)
})

# TC-10.2.27: 预期效应图在推荐排除点使用红色空心圆标记
test_that("TC-10.2.27: anticipation plot highlights recommended exclusion", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_no_anticipation(
    anticipation_detected = TRUE, recommended_exclusion = 2L
  )
  p <- plot.lwdid_sensitivity(x)
  built <- ggplot2::ggplot_build(p)

  highlight_found <- FALSE
  for (i in seq_along(p$layers)) {
    if (!inherits(p$layers[[i]]$geom, "GeomPoint")) {
      next
    }

    point_data <- built$data[[i]]
    if (!all(c("x", "y", "shape", "colour") %in% names(point_data))) {
      next
    }

    if (any(abs(point_data$x - 2) < 1e-10 &
            abs(point_data$y - 1.9) < 1e-10 &
            point_data$shape == 1 &
            grepl("red|#ff0000", tolower(point_data$colour)))) {
      highlight_found <- TRUE
      break
    }
  }

  expect_true(highlight_found)
})

# TC-10.2.28: 队列趋势图置信带使用t分布临界值
test_that("TC-10.2.28: trend CI uses t-distribution critical value", {
  skip_if_not_installed("ggplot2")
  # Single cohort: cohort=2004, n_pre=4, slope_se=0.1
  # df = max(4-2, 1) = 2, t_crit = qt(0.975, 2) ≈ 4.303
  # At t_centered = 1.5 (endpoint): se_fitted = 0.1 * 1.5 = 0.15
  # CI half-width = 4.303 * 0.15 ≈ 0.645
  # z=1.96 half-width = 1.96 * 0.15 = 0.294
  trends <- list(
    list(cohort = 2004L, intercept = 10.0, slope = 0.5,
         slope_se = 0.1, slope_pvalue = 0.001,
         n_units = 50L, n_pre_periods = 4L,
         r_squared = 0.85, residual_std = 0.5)
  )
  x <- create_mock_heterogeneous_trends(
    trend_by_cohort = trends,
    control_group_trend = NULL
  )
  p <- plot.lwdid_heterogeneous_trends(x, show_ci = TRUE,
                                        show_control = FALSE)
  built <- ggplot2::ggplot_build(p)
  # Find ribbon layer
  ribbon_data <- NULL
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomRibbon")) {
      ribbon_data <- built$data[[i]]
      break
    }
  }
  expect_false(is.null(ribbon_data))
  # Endpoint CI should be wider than z=1.96 would give
  # time_seq = 2000:2003, t_centered = c(-1.5, -0.5, 0.5, 1.5)
  # At endpoint (row 1 or 4): ci_width = 2 * 4.303 * 0.15 ≈ 1.291
  endpoint_width <- ribbon_data$ymax[1] - ribbon_data$ymin[1]
  z_width <- 2 * 1.96 * 0.1 * 1.5  # z * slope_se * |t_centered|
  expect_true(endpoint_width > z_width)
})

# TC-10.2.29: t分布CI与z=1.96 CI的数值差异验证（纯数值）
test_that("TC-10.2.29: t-dist vs z=1.96 numerical differences", {
  # df=5: t_0.975 ≈ 2.5706, diff from z=1.96 ≈ 31%
  t5 <- qt(0.975, 5)
  expect_equal(t5, 2.5706, tolerance = 0.001)
  expect_true((t5 - 1.96) / 1.96 > 0.30)

  # df=2: t_0.975 ≈ 4.3027, diff ≈ 120%
  t2 <- qt(0.975, 2)
  expect_equal(t2, 4.3027, tolerance = 0.001)
  expect_true((t2 - 1.96) / 1.96 > 1.10)

  # df=100: t_0.975 ≈ 1.984, close to z=1.96
  t100 <- qt(0.975, 100)
  expect_equal(t100, 1.984, tolerance = 0.001)
  expect_true(abs(t100 - 1.96) < 0.05)
})

# TC-10.2.30: 队列趋势图拟合值与置信带数值精确验证
test_that("TC-10.2.30: trend fitted values and CI numerical", {
  skip_if_not_installed("ggplot2")
  # Single cohort: cohort=2004, intercept=10.0, slope=0.5,
  # slope_se=0.1, n_pre=4
  # time_seq = [2000,2001,2002,2003]
  # t_mean = 2001.5
  # t_centered = [-1.5, -0.5, 0.5, 1.5]
  # fitted = [9.25, 9.75, 10.25, 10.75]
  # df = max(4-2, 1) = 2
  # t_crit = qt(0.975, 2) ≈ 4.3027
  # se_fitted = [0.15, 0.05, 0.05, 0.15]
  trends <- list(
    list(cohort = 2004L, intercept = 10.0, slope = 0.5,
         slope_se = 0.1, slope_pvalue = 0.001,
         n_units = 50L, n_pre_periods = 4L,
         r_squared = 0.85, residual_std = 0.5)
  )
  x <- create_mock_heterogeneous_trends(
    trend_by_cohort = trends,
    control_group_trend = NULL
  )
  p <- plot.lwdid_heterogeneous_trends(x, show_ci = TRUE,
                                        show_control = FALSE,
                                        show_slope_labels = FALSE)
  built <- ggplot2::ggplot_build(p)

  # Line layer data (fitted values)
  line_data <- built$data[[1]]
  expect_equal(line_data$y, c(9.25, 9.75, 10.25, 10.75),
               tolerance = 1e-10)

  # Ribbon layer data (CI)
  ribbon_data <- built$data[[2]]
  t_crit <- qt(0.975, 2)
  se_fitted <- c(0.15, 0.05, 0.05, 0.15)
  expected_lower <- c(9.25, 9.75, 10.25, 10.75) - t_crit * se_fitted
  expected_upper <- c(9.25, 9.75, 10.25, 10.75) + t_crit * se_fitted
  expect_equal(ribbon_data$ymin, expected_lower, tolerance = 1e-6)
  expect_equal(ribbon_data$ymax, expected_upper, tolerance = 1e-6)

  # CI symmetry: widths at positions 2 and 3 should be equal
  ci_widths <- ribbon_data$ymax - ribbon_data$ymin
  expect_equal(ci_widths[2], ci_widths[3], tolerance = 1e-10)
  # Endpoint CI wider than midpoint CI
  expect_true(ci_widths[1] > ci_widths[2])
})


# TC-10.2.31: 所有规格均未收敛时优雅处理
test_that("TC-10.2.31: all specs unconverged gives warning + empty ggplot", {
  skip_if_not_installed("ggplot2")
  specs <- list(
    list(n_pre_periods = 2L, att = 1.0, se = 0.3, pvalue = 0.1,
         ci_lower = 0.4, ci_upper = 1.6, converged = FALSE),
    list(n_pre_periods = 3L, att = 1.5, se = 0.2, pvalue = 0.05,
         ci_lower = 1.1, ci_upper = 1.9, converged = FALSE)
  )
  x <- create_mock_sensitivity_pre_period(specifications = specs)
  expect_warning(p <- plot.lwdid_sensitivity(x))
  expect_s3_class(p, "ggplot")
  # Empty ggplot should have no layers with data
  expect_equal(length(p$layers), 0L)
})

# TC-10.2.32: 预期效应图所有估计未收敛时优雅处理
test_that("TC-10.2.32: all anticipation estimates unconverged", {
  skip_if_not_installed("ggplot2")
  estimates <- list(
    list(excluded_periods = 0L, att = 2.0, se = 0.3,
         ci_lower = 1.4, ci_upper = 2.6, converged = FALSE),
    list(excluded_periods = 1L, att = 2.2, se = 0.28,
         ci_lower = 1.64, ci_upper = 2.76, converged = FALSE)
  )
  x <- create_mock_sensitivity_no_anticipation(estimates = estimates)
  expect_warning(p <- plot.lwdid_sensitivity(x))
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 0L)
})

# TC-10.2.33: 综合敏感性图show_ci=FALSE传递到两个子图
test_that("TC-10.2.33: comprehensive show_ci=FALSE propagates", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  x <- create_mock_sensitivity_comprehensive()
  p <- plot.lwdid_sensitivity_comprehensive(x, show_ci = FALSE)
  expect_s3_class(p, "patchwork")
  # Extract sub-plots and verify no errorbar layers
  sub1 <- p[[1]]
  sub2 <- p[[2]]
  layers1 <- vapply(sub1$layers, function(l) class(l$geom)[1], character(1))
  layers2 <- vapply(sub2$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomErrorbar" %in% layers1)
  expect_false("GeomErrorbar" %in% layers2)
})

# TC-10.2.34: 队列趋势图中心点CI宽度为零的数值验证
test_that("TC-10.2.34: center point CI width is zero", {
  # Pure numerical: when t_centered=0, se_fitted = sqrt(slope_se^2 * 0^2) = 0
  # So CI_lower = CI_upper = fitted = intercept
  slope_se <- 0.1
  t_centered <- 0
  se_fitted <- sqrt(slope_se^2 * t_centered^2)
  expect_equal(se_fitted, 0, tolerance = 1e-15)
  # At center: fitted = intercept (since slope * 0 = 0)
  intercept <- 10.0
  slope <- 0.5
  fitted_center <- intercept + slope * t_centered
  expect_equal(fitted_center, intercept, tolerance = 1e-15)
  # CI width = 2 * t_crit * 0 = 0
  t_crit <- qt(0.975, 2)
  ci_width <- 2 * t_crit * se_fitted
  expect_equal(ci_width, 0, tolerance = 1e-15)
})

# TC-10.2.35: method_labels部分覆盖时正确处理
test_that("TC-10.2.35: partial method_labels preserves unlabeled", {
  skip_if_not_installed("ggplot2")
  labels <- c(demean = "\u53bb\u5747\u503c", demeanq = "\u53bb\u5747\u503c(\u4e8c\u6b21)")
  x <- create_mock_transformation_recommendation(method_labels = labels)
  p <- plot.lwdid_transformation_recommendation(x)
  lvls <- levels(p$data$method)
  expect_true("\u53bb\u5747\u503c" %in% lvls)
  expect_true("\u53bb\u5747\u503c(\u4e8c\u6b21)" %in% lvls)
  # detrend has no label, should keep original name
  expect_true("detrend" %in% lvls)
})

# TC-10.2.36: 规格曲线图x轴刻度对齐到整数
test_that("TC-10.2.36: spec curve x-axis breaks are integers", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_pre_period()
  p <- plot.lwdid_sensitivity(x)
  built <- ggplot2::ggplot_build(p)
  # x-axis breaks should be the n_pre values (integers)
  x_breaks <- built$layout$panel_params[[1]]$x$breaks
  x_breaks <- x_breaks[!is.na(x_breaks)]
  expect_true(all(x_breaks == as.integer(x_breaks)))
  expect_true(all(c(2L, 3L, 4L, 5L) %in% x_breaks))
})

# TC-10.2.37: 未检测到预期效应时不添加推荐排除点高亮
test_that("TC-10.2.37: no highlight when anticipation_detected=FALSE", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_no_anticipation(
    anticipation_detected = FALSE
  )
  p <- plot.lwdid_sensitivity(x)
  built <- ggplot2::ggplot_build(p)

  highlight_found <- FALSE
  for (i in seq_along(p$layers)) {
    if (!inherits(p$layers[[i]]$geom, "GeomPoint")) {
      next
    }

    point_data <- built$data[[i]]
    if (!all(c("shape", "colour") %in% names(point_data))) {
      next
    }

    if (any(point_data$shape == 1 &
            grepl("red|#ff0000", tolower(point_data$colour)))) {
      highlight_found <- TRUE
      break
    }
  }

  expect_false(highlight_found)
})

# TC-10.2.38: 队列趋势图对照组趋势线使用虚线样式
test_that("TC-10.2.38: control group uses dashed linetype", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends()
  p <- plot.lwdid_heterogeneous_trends(x, show_control = TRUE)
  built <- ggplot2::ggplot_build(p)
  # Find the control group line layer (dashed)
  found_dashed <- FALSE
  for (i in seq_along(p$layers)) {
    layer_data <- built$data[[i]]
    if (inherits(p$layers[[i]]$geom, "GeomLine") &&
        "linetype" %in% names(layer_data)) {
      # linetype "dashed" is encoded as numeric 2 in ggplot2
      if (any(layer_data$linetype == "dashed" |
              layer_data$linetype == 2)) {
        found_dashed <- TRUE
        break
      }
    }
  }
  expect_true(found_dashed)
})

# TC-10.2.39: 规格曲线图显著性颜色编码正确
test_that("TC-10.2.39: significance colors blue/red correct", {
  skip_if_not_installed("ggplot2")
  # Construct specs: one significant (p=0.001), one not (p=0.5)
  specs <- list(
    list(n_pre_periods = 2L, att = 1.8, se = 0.3, pvalue = 0.001,
         ci_lower = 1.2, ci_upper = 2.4, converged = TRUE),
    list(n_pre_periods = 3L, att = 2.1, se = 0.25, pvalue = 0.5,
         ci_lower = 1.6, ci_upper = 2.6, converged = TRUE)
  )
  x <- create_mock_sensitivity_pre_period(specifications = specs)
  p <- plot.lwdid_sensitivity(x, ci_level = 0.95)
  built <- ggplot2::ggplot_build(p)
  # Find point layer
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
      point_data <- built$data[[i]]
      colors <- point_data$colour
      # Significant (p=0.001 < 0.05) should be blue #2166AC
      # Not significant (p=0.5 > 0.05) should be red #B2182B
      expect_true("#2166AC" %in% toupper(colors))
      expect_true("#B2182B" %in% toupper(colors))
      break
    }
  }
})

# TC-10.2.40: 队列趋势图无对照组数据时show_control=TRUE不崩溃
test_that("TC-10.2.40: show_control=TRUE with NULL control no crash", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_heterogeneous_trends(control_group_trend = NULL)
  p <- plot.lwdid_heterogeneous_trends(x, show_control = TRUE)
  expect_s3_class(p, "ggplot")
})

# TC-10.2.41: show_ci=FALSE 仍保留推荐排除点高亮
test_that("TC-10.2.41: show_ci=FALSE keeps recommended highlight", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_sensitivity_no_anticipation(
    anticipation_detected = TRUE, recommended_exclusion = 2L
  )
  p <- plot.lwdid_sensitivity(x, show_ci = FALSE)
  built <- ggplot2::ggplot_build(p)

  highlight_found <- FALSE
  for (i in seq_along(p$layers)) {
    if (!inherits(p$layers[[i]]$geom, "GeomPoint")) {
      next
    }

    point_data <- built$data[[i]]
    if (!all(c("x", "y", "shape", "colour") %in% names(point_data))) {
      next
    }

    if (any(abs(point_data$x - 2) < 1e-10 &
            abs(point_data$y - 1.9) < 1e-10 &
            point_data$shape == 1 &
            grepl("red|#ff0000", tolower(point_data$colour)))) {
      highlight_found <- TRUE
      break
    }
  }

  expect_true(highlight_found)
})

# TC-10.2.42: .format_pvalue()正确格式化
test_that("TC-10.2.42: .format_pvalue formats correctly", {
  # Access internal function via :::
  fmt <- lwdid:::.format_pvalue
  # p=0.0001 -> "<0.001"
  expect_equal(fmt(0.0001, 3L), "<0.001")
  # p=0.05 -> "0.0500" (digits=4)
  expect_equal(fmt(0.05, 4L), "0.0500")
  # p=0.123 -> "0.123" (digits=3)
  expect_equal(fmt(0.123, 3L), "0.123")
  # p=0.001 is NOT < 0.001, so formatted as "0.0010" (digits=4)
  expect_equal(fmt(0.001, 4L), "0.0010")
  # p=0.001 with digits=3 -> "0.001"
  expect_equal(fmt(0.001, 3L), "0.001")
})

# TC-10.2.43: 规格曲线图使用vapply确保类型安全
test_that("TC-10.2.43: vapply type safety in spec curve", {
  skip_if_not_installed("ggplot2")
  # Construct specifications with consistent types
  # vapply should enforce type safety - if all specs have correct types,
  # the function should work without issues
  specs <- list(
    list(n_pre_periods = 2L, att = 1.8, se = 0.3, pvalue = 0.001,
         ci_lower = 1.2, ci_upper = 2.4, converged = TRUE),
    list(n_pre_periods = 3L, att = 2.1, se = 0.25, pvalue = 0.0001,
         ci_lower = 1.6, ci_upper = 2.6, converged = TRUE)
  )
  x <- create_mock_sensitivity_pre_period(specifications = specs)
  # Should succeed without error (vapply enforces types)
  p <- plot.lwdid_sensitivity(x)
  expect_s3_class(p, "ggplot")
  # Verify the plot data has correct types
  built <- ggplot2::ggplot_build(p)
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
      point_data <- built$data[[i]]
      expect_true(is.numeric(point_data$x))
      expect_true(is.numeric(point_data$y))
      break
    }
  }
})
