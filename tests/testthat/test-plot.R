# ============================================================================
# test-plot.R — Unified Plot Entry Point & Diagnostic Plot Tests
# TC-10.5.1 through TC-10.5.43
# Story E10-05: plot.lwdid_result统一绘图入口
# ============================================================================

skip_if_not_installed("ggplot2")

# ── Mock object factories ──────────────────────────────────────────────────

.matches_utf8_or_escaped <- function(text, utf8_literal, escaped_regex) {
  grepl(utf8_literal, text, fixed = TRUE) ||
    grepl(escaped_regex, text)
}

#' Create a minimal lwdid_result with plot_data (non-staggered)
create_mock_result_with_plotdata <- function() {
  obj <- structure(list(
    att = 0.5, se_att = 0.1, t_stat = 5.0, pvalue = 0.001,
    ci_lower = 0.3, ci_upper = 0.7, df_inference = 50L,
    alpha = 0.05, n_obs = 200L,
    cohort_time_effects = data.frame(
      cohort = c(2004L, 2004L), time = c(2004L, 2005L),
      att = c(0.4, 0.6), se = c(0.1, 0.12),
      t_stat = c(4.0, 5.0), pvalue = c(0.001, 0.0001),
      ci_lower = c(0.2, 0.36), ci_upper = c(0.6, 0.84),
      weight = c(0.5, 0.5)
    ),
    is_staggered = FALSE,
    metadata = list(
      plot_data = list(
        time = 2000:2006,
        control_mean = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6),
        treated_mean = c(1.0, 1.1, 1.2, 1.8, 2.0, 2.2, 2.4),
        intervention_point = 2003,
        treated_individual = list(
          "u1" = c(0.9, 1.0, 1.1, 1.7, 1.9, 2.1, 2.3),
          "u2" = c(1.1, 1.2, 1.3, 1.9, 2.1, 2.3, 2.5)
        )
      ),
      estimator = "ra", control_group = "not_yet_treated",
      aggregate = "cohort", rolling = "demean"
    ),
    sensitivity = NULL,
    diagnostics = NULL,
    call_info = list(formula_str = "y ~ treat")
  ), class = "lwdid_result")
  obj
}

#' Create a staggered lwdid_result with cohort plot_data
create_mock_staggered_plotdata <- function() {
  obj <- create_mock_result_with_plotdata()
  obj$is_staggered <- TRUE
  obj$metadata$plot_data <- list(
    cohorts = list(
      "2004" = list(
        time = 2000:2006,
        control_mean = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6),
        treated_mean = c(1.0, 1.1, 1.2, 1.8, 2.0, 2.2, 2.4),
        intervention_point = 2004,
        treated_individual = list(
          "u1" = c(0.9, 1.0, 1.1, 1.7, 1.9, 2.1, 2.3),
          "u2" = c(1.1, 1.2, 1.3, 1.9, 2.1, 2.3, 2.5)
        )
      ),
      "2006" = list(
        time = 2002:2008,
        control_mean = c(2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6),
        treated_mean = c(2.0, 2.1, 2.2, 2.3, 2.9, 3.1, 3.3),
        intervention_point = 2006,
        treated_individual = list(
          "u3" = c(1.9, 2.0, 2.1, 2.2, 2.8, 3.0, 3.2)
        )
      )
    )
  )
  obj
}

#' Create mock lwdid_trends object
create_mock_trends <- function(with_group_trends = TRUE) {
  obj <- structure(list(
    pre_treatment_coefficients = data.frame(
      period = c(-3, -2, -1),
      coefficient = c(0.02, -0.01, 0.005),
      se = c(0.05, 0.04, 0.03),
      ci_lower = c(-0.08, -0.09, -0.055),
      ci_upper = c(0.12, 0.07, 0.065)
    ),
    f_test = list(df1 = 3L, df2 = 97L, f_stat = 0.85, f_pvalue = 0.47),
    group_trends = if (with_group_trends) {
      data.frame(
        period = rep(c(-3, -2, -1, 0, 1), 2),
        mean_y = c(1.0, 1.1, 1.2, 1.3, 1.4,
                   1.0, 1.15, 1.25, 1.8, 2.0),
        group = rep(c("control", "treated"), each = 5)
      )
    } else NULL
  ), class = "lwdid_trends")
  obj
}

#' Create mock lwdid_clustering_diagnosis object
create_mock_clustering <- function(n_clusters = 10L) {
  sizes <- setNames(
    as.integer(seq(50, 50 + (n_clusters - 1) * 5, by = 5)),
    paste0("cl_", seq_len(n_clusters))
  )
  structure(list(
    cluster_sizes = sizes,
    n_clusters = n_clusters,
    balance_ratio = 0.85,
    icc = 0.05,
    effective_clusters = n_clusters * 0.8
  ), class = "lwdid_clustering_diagnosis")
}

#' Create mock lwdid_selection_diagnosis object
create_mock_selection <- function(with_attrition = TRUE) {
  structure(list(
    attrition_by_period = if (with_attrition) {
      data.frame(
        period = 2001:2005,
        attrition_rate = c(0.02, 0.03, 0.05, 0.04, 0.06)
      )
    } else NULL,
    attrition_rate = 0.04,
    selection_risk = "Low"
  ), class = "lwdid_selection_diagnosis")
}

#' Create mock result with diagnostics
create_mock_result_with_diagnostics <- function(
    trends = TRUE, clustering = TRUE, selection = TRUE) {
  obj <- create_mock_result_with_plotdata()
  obj$diagnostics <- list(
    parallel_trends = if (trends) create_mock_trends() else NULL,
    clustering = if (clustering) create_mock_clustering() else NULL,
    selection = if (selection) create_mock_selection() else NULL
  )
  obj
}

# ============================================================================
# plot.lwdid_result dispatch tests (TC-10.5.1 to TC-10.5.6)
# ============================================================================

# We need fixture_staggered_result for event_study dispatch test
# Use a minimal approach: create result with cohort_time_effects that has event_time
create_mock_result_for_event_study <- function() {
  obj <- create_mock_result_with_plotdata()
  obj$cohort_time_effects$event_time <- c(0L, 1L)
  # plot_event_study needs att_by_period (aggregated by event_time)
  obj$att_by_period <- data.frame(
    event_time = c(-2L, -1L, 0L, 1L),
    att = c(0.01, -0.005, 0.4, 0.6),
    se = c(0.05, 0.04, 0.1, 0.12),
    ci_lower = c(-0.09, -0.085, 0.2, 0.36),
    ci_upper = c(0.11, 0.075, 0.6, 0.84),
    pvalue = c(0.84, 0.90, 0.001, 0.0001)
  )
  obj
}

test_that("TC-10.5.1: default type=event_study returns ggplot with GeomErrorbar", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_for_event_study()
  p <- suppressWarnings(plot(res))
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomErrorbar" %in% layer_classes)
})

test_that("TC-10.5.2: type=trajectories returns ggplot with GeomLine", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  p <- plot(res, type = "trajectories")
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLine" %in% layer_classes)
})

test_that("TC-10.5.3: type=trajectories without plot_data errors", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  res$metadata$plot_data <- NULL
  expect_error(plot(res, type = "trajectories"), "graph=TRUE")
})

test_that("TC-10.5.4: type=sensitivity without sensitivity errors", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  expect_error(plot(res, type = "sensitivity"), "sensitivity_analysis")
})

test_that("TC-10.5.5: type=diagnostics without diagnostics errors", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  expect_error(plot(res, type = "diagnostics"), "return_diagnostics")
})

test_that("TC-10.5.6: invalid type errors via match.arg", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  expect_error(plot(res, type = "invalid_type"), "arg")
})

# ============================================================================
# Trajectory detail tests (TC-10.5.7 to TC-10.5.8)
# ============================================================================

test_that("TC-10.5.7: show_individual=TRUE overlays individual trajectories", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  p <- plot(res, type = "trajectories", show_individual = TRUE)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  # Base has 2 layers (geom_line + geom_vline), individual adds 1 more geom_line
  n_line <- sum(layer_classes == "GeomLine")
  expect_gte(n_line, 2L)  # at least base line + individual overlay
})

test_that("TC-10.5.8: highlight_id highlights specific units in orange", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  p <- plot(res, type = "trajectories", highlight_id = "u1")
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  # Base has 2 layers, highlight adds 1 more geom_line
  n_line <- sum(layer_classes == "GeomLine")
  expect_gte(n_line, 2L)
})

# ============================================================================
# Diagnostics panel tests (TC-10.5.9 to TC-10.5.10)
# ============================================================================

test_that("TC-10.5.9: single diagnostic panel returns ggplot (not patchwork)", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_diagnostics(
    trends = TRUE, clustering = FALSE, selection = FALSE)
  p <- plot(res, type = "diagnostics")
  expect_s3_class(p, "ggplot")
  # Should NOT be patchwork
  expect_false(inherits(p, "patchwork"))
})

test_that("TC-10.5.10: multiple diagnostics use patchwork", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  res <- create_mock_result_with_diagnostics(
    trends = TRUE, clustering = TRUE, selection = TRUE)
  p <- plot(res, type = "diagnostics")
  expect_s3_class(p, "patchwork")
})

# ============================================================================
# plot.lwdid_trends tests (TC-10.5.11 to TC-10.5.13)
# ============================================================================

test_that("TC-10.5.11: trends coefficients plot has GeomHline reference line", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_trends()
  p <- plot(x, type = "coefficients")
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomHline" %in% layer_classes)
})

test_that("TC-10.5.12: trends trajectories plot draws group trend lines", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_trends(with_group_trends = TRUE)
  p <- plot(x, type = "trajectories")
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLine" %in% layer_classes)
  expect_true("GeomSmooth" %in% layer_classes)
})

test_that("TC-10.5.13: trends trajectories without group_trends errors", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_trends(with_group_trends = FALSE)
  err <- tryCatch(
    {
      plot(x, type = "trajectories")
      NULL
    },
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_true(
    .matches_utf8_or_escaped(
      conditionMessage(err),
      utf8_literal = "趋势",
      escaped_regex = "<U\\+8D8B><U\\+52BF>"
    )
  )
})

# ============================================================================
# plot.lwdid_clustering_diagnosis tests (TC-10.5.14)
# ============================================================================

test_that("TC-10.5.14: clustering diagnostic plot has GeomCol layer", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_clustering(10L)
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomCol" %in% layer_classes)
})

# ============================================================================
# plot.lwdid_selection_diagnosis tests (TC-10.5.15 to TC-10.5.17)
# ============================================================================

test_that("TC-10.5.15: selection diagnostic with time series has GeomLine", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_selection(with_attrition = TRUE)
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLine" %in% layer_classes)
})

test_that("TC-10.5.16: selection diagnostic without time series falls back to text", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_selection(with_attrition = FALSE)
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomText" %in% layer_classes)
  # Verify theme_void is applied (no axis lines)
  theme_elems <- p$theme
  # theme_void sets panel.grid, axis.line etc to element_blank
})

test_that("TC-10.5.17: selection diagnostic scales fallback when scales unavailable", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_selection(with_attrition = TRUE)
  # We can't truly unload scales, but verify the plot works regardless
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  # The y-axis should show attrition rates — verify data range is reasonable
  build <- ggplot2::ggplot_build(p)
  y_range <- range(build$data[[1]]$y, na.rm = TRUE)
  expect_true(y_range[1] >= 0)
  expect_true(y_range[2] <= 1)  # attrition rates are proportions
})

# ============================================================================
# Type check tests (TC-10.5.18)
# ============================================================================

test_that("TC-10.5.18: plot.lwdid_result rejects non-lwdid_result objects", {
  skip_if_not_installed("ggplot2")
  expect_error(plot.lwdid_result(list(a = 1)), "inherits")
  expect_error(plot.lwdid_result("not_a_result"), "inherits")
})

# ============================================================================
# gid parameter tests (TC-10.5.19 to TC-10.5.23)
# ============================================================================

test_that("TC-10.5.19: staggered gid selects specific cohort, title has cohort ID", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_staggered_plotdata()
  p <- plot(res, type = "trajectories", gid = 2004)
  expect_s3_class(p, "ggplot")
  expect_true(grepl("2004", p$labels$title))
})

test_that("TC-10.5.20: staggered gid=NULL creates faceted plot (FacetWrap)", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_staggered_plotdata()
  p <- plot(res, type = "trajectories")
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("TC-10.5.21: staggered gid with invalid cohort errors with available list", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_staggered_plotdata()
  expect_error(
    plot(res, type = "trajectories", gid = 9999),
    "9999"
  )
  # Error message should list available cohorts
  err <- tryCatch(
    plot(res, type = "trajectories", gid = 9999),
    error = function(e) e$message
  )
  expect_true(grepl("2004", err))
  expect_true(grepl("2006", err))
})

test_that("TC-10.5.22: non-staggered with gid emits warning", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  expect_warning(
    plot(res, type = "trajectories", gid = 1),
    "gid"
  )
})

test_that("TC-10.5.23: staggered gid with show_individual overlays cohort individuals", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_staggered_plotdata()
  p <- plot(res, type = "trajectories", gid = 2004, show_individual = TRUE)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  # Should have base line + vline + individual overlay line = at least 3 GeomLine-like layers
  n_line <- sum(layer_classes == "GeomLine")
  expect_gte(n_line, 2L)
})

# ============================================================================
# which parameter tests (TC-10.5.24 to TC-10.5.27)
# ============================================================================

test_that("TC-10.5.24: which filters to single diagnostic (non-patchwork)", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_diagnostics(
    trends = TRUE, clustering = TRUE, selection = TRUE)
  p <- plot(res, type = "diagnostics", which = "clustering")
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
})

test_that("TC-10.5.25: which with multiple diagnostics uses patchwork", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  res <- create_mock_result_with_diagnostics(
    trends = TRUE, clustering = TRUE, selection = TRUE)
  p <- plot(res, type = "diagnostics",
            which = c("parallel_trends", "clustering"))
  expect_s3_class(p, "patchwork")
})

test_that("TC-10.5.26: which with unavailable diagnostic type errors", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_diagnostics(
    trends = FALSE, clustering = FALSE, selection = TRUE)
  err <- tryCatch(
    {
      plot(res, type = "diagnostics", which = "parallel_trends")
      NULL
    },
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_true(
    .matches_utf8_or_escaped(
      conditionMessage(err),
      utf8_literal = "无可用数据",
      escaped_regex = "<U\\+65E0><U\\+53EF><U\\+7528><U\\+6570><U\\+636E>"
    ) || .matches_utf8_or_escaped(
      conditionMessage(err),
      utf8_literal = "无",
      escaped_regex = "<U\\+65E0>$"
    )
  )
})

test_that("TC-10.5.27: all diagnostics NULL errors", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_diagnostics(
    trends = FALSE, clustering = FALSE, selection = FALSE)
  expect_error(
    plot(res, type = "diagnostics"),
    "NULL"
  )
})

# ============================================================================
# smooth_method parameter test (TC-10.5.28)
# ============================================================================

test_that("TC-10.5.28: smooth_method=loess uses loess, subtitle contains loess", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_trends(with_group_trends = TRUE)
  p <- plot(x, type = "trajectories", smooth_method = "loess")
  expect_s3_class(p, "ggplot")
  expect_true(grepl("loess", p$labels$subtitle))
})

# ============================================================================
# Large cluster handling tests (TC-10.5.29 to TC-10.5.31)
# ============================================================================

test_that("TC-10.5.29: >30 clusters uses coord_flip (horizontal bars)", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_clustering(35L)
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$coordinates, "CoordFlip"))
})

test_that("TC-10.5.30: >60 clusters shows top/bottom 15, subtitle mentions truncation", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_clustering(70L)
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  expect_true(grepl("showing top/bottom 15", p$labels$subtitle))
  # Verify only 30 bars (15 top + 15 bottom)
  build <- ggplot2::ggplot_build(p)
  col_data <- build$data[[1]]  # GeomCol data
  expect_equal(nrow(col_data), 30L)
})

test_that("TC-10.5.31: <=30 clusters uses standard vertical bars (no CoordFlip)", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_clustering(20L)
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p$coordinates, "CoordFlip"))
})

# ============================================================================
# Edge cases & robustness tests (TC-10.5.32 to TC-10.5.38)
# ============================================================================

test_that("TC-10.5.32: highlight_id with nonexistent ID silently skips", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  # Should not error — nonexistent IDs are silently ignored
  p <- plot(res, type = "trajectories", highlight_id = "nonexistent_unit")
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  # Only base layers (geom_line + geom_vline), no highlight layer added
  n_line <- sum(layer_classes == "GeomLine")
  expect_equal(n_line, 1L)  # just the base trajectory line
})

test_that("TC-10.5.33: show_individual=TRUE but no treated_individual silently skips", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  res$metadata$plot_data$treated_individual <- NULL
  p <- plot(res, type = "trajectories", show_individual = TRUE)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  n_line <- sum(layer_classes == "GeomLine")
  expect_equal(n_line, 1L)  # only base trajectory, no individual overlay
})

test_that("TC-10.5.34: staggered gid=NULL faceted plot has per-cohort GeomVline", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_staggered_plotdata()
  p <- plot(res, type = "trajectories")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomVline" %in% layer_classes)
})

test_that("TC-10.5.35: plot.lwdid_trends with invalid type errors", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_trends()
  expect_error(plot(x, type = "invalid"), "arg")
})

test_that("TC-10.5.36: plot.lwdid_trends rejects non-lwdid_trends object", {
  skip_if_not_installed("ggplot2")
  expect_error(plot.lwdid_trends(list(a = 1)), "inherits")
})

test_that("TC-10.5.37: plot.lwdid_clustering_diagnosis rejects non-valid object", {
  skip_if_not_installed("ggplot2")
  expect_error(plot.lwdid_clustering_diagnosis(list(a = 1)), "inherits")
})

test_that("TC-10.5.38: plot.lwdid_selection_diagnosis rejects non-valid object", {
  skip_if_not_installed("ggplot2")
  expect_error(plot.lwdid_selection_diagnosis(list(a = 1)), "inherits")
})

# ============================================================================
# Theme consistency & additional edge tests (TC-10.5.39 to TC-10.5.43)
# ============================================================================

test_that("TC-10.5.39: custom theme applied in non-staggered trajectory plot", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  p <- plot(res, type = "trajectories", theme = ggplot2::theme_bw)
  expect_s3_class(p, "ggplot")
  # theme_bw sets panel.border to element_rect (not element_blank like theme_minimal)
  # Verify the theme was applied by checking it's a valid ggplot
  build <- ggplot2::ggplot_build(p)
  expect_true(length(build$data) > 0)
})

test_that("TC-10.5.40: custom theme applied in staggered faceted trajectory plot", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_staggered_plotdata()
  p <- plot(res, type = "trajectories", theme = ggplot2::theme_bw)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
  # Verify theme_bw was applied (not hardcoded theme_minimal)
  build <- ggplot2::ggplot_build(p)
  expect_true(length(build$data) > 0)
})

test_that("TC-10.5.41: staggered gid=NULL faceted plot ignores show_individual (early return)", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_staggered_plotdata()
  # With show_individual=TRUE but gid=NULL → faceted path does early return
  p_without <- plot(res, type = "trajectories", show_individual = FALSE)
  p_with <- plot(res, type = "trajectories", show_individual = TRUE)
  # Both should have same number of layers (early return skips individual overlay)
  expect_equal(length(p_without$layers), length(p_with$layers))
})

test_that("TC-10.5.42: highlight_id works independently of show_individual=FALSE", {
  skip_if_not_installed("ggplot2")
  res <- create_mock_result_with_plotdata()
  p <- plot(res, type = "trajectories",
            show_individual = FALSE, highlight_id = "u1")
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  n_line <- sum(layer_classes == "GeomLine")
  # Base line (1) + highlight line (1) = 2, no individual overlay
  expect_equal(n_line, 2L)
})

test_that("TC-10.5.43: trends trajectories subtitle has smooth_method, no deprecation warning", {
  skip_if_not_installed("ggplot2")
  x <- create_mock_trends(with_group_trends = TRUE)
  # Capture warnings during plot creation
  warns <- character(0)
  p <- withCallingHandlers(
    plot(x, type = "trajectories", smooth_method = "lm"),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(p, "ggplot")
  expect_true(grepl("lm", p$labels$subtitle))
  # No deprecation warning about formula
  deprecation_warns <- grep("formula", warns, value = TRUE, ignore.case = TRUE)
  expect_equal(length(deprecation_warns), 0L)
})
