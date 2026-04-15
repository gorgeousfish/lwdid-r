# ===========================================================================
# plot_diagnostics.R — Sensitivity analysis and trend diagnostic plots
#
# Implements specification curve plots, anticipated effect sensitivity plots,
# cohort trend plots, and transformation recommendation score plots.
# All plots built with ggplot2 layer grammar; composite plots use patchwork.
# ===========================================================================

#' @title Sensitivity and Diagnostic Plots
#' @description S3 plot methods for sensitivity analysis and diagnostic objects.
#' @name plot_diagnostics
#' @family lwdid-visualization
NULL


# NOTE: .format_pvalue() is defined in results.R (Story 10.4)
# and shared across all modules. Do not duplicate here.


# ---- Internal: specification curve plot ----

.plot_specification_curve <- function(x, ci_level = 0.95, show_threshold = TRUE,
                                      show_ci = TRUE, ...) {
  # --- 1. Data preparation ---
  alpha <- 1 - ci_level

  # Filter converged specifications
  converged_mask <- vapply(x$specifications, function(s) s$converged, logical(1))
  specs <- x$specifications[converged_mask]

  # Empty data guard
  if (length(specs) == 0L) {
    warning("\u6240\u6709\u89c4\u683c\u5747\u672a\u6536\u655b\uff0c\u65e0\u6cd5\u7ed8\u5236\u6709\u6548\u56fe\u5f62")
    p <- ggplot2::ggplot() +
      ggplot2::labs(
        title = "\u524d\u671f\u7a33\u5065\u6027\u89c4\u683c\u66f2\u7ebf",
        subtitle = "\u65e0\u6536\u655b\u89c4\u683c"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_text(color = "red")
      )
    return(p)
  }

  # Build plot_data data.frame using vapply
  plot_data <- data.frame(
    n_pre       = vapply(specs, function(s) s$n_pre_periods, integer(1)),
    att         = vapply(specs, function(s) s$att, numeric(1)),
    ci_lower    = vapply(specs, function(s) s$ci_lower, numeric(1)),
    ci_upper    = vapply(specs, function(s) s$ci_upper, numeric(1)),
    significant = vapply(specs, function(s) s$pvalue < alpha, logical(1))
  )

  # --- 2. Build ggplot ---
  baseline_att <- x$baseline_spec$att
  threshold    <- x$robustness_threshold

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = n_pre, y = att))

  # Threshold band
  if (show_threshold) {
    p <- p + ggplot2::annotate(
      "rect",
      xmin  = -Inf,
      xmax  = Inf,
      ymin  = baseline_att - threshold * abs(baseline_att),
      ymax  = baseline_att + threshold * abs(baseline_att),
      fill  = "gray90",
      alpha = 0.5
    )
  }

  # Reference lines
  p <- p +
    ggplot2::geom_hline(yintercept = baseline_att, linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "gray80")

  # Confidence intervals
  if (show_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.3,
      color = "gray40"
    )
  }

  # Point estimates with significance coloring
  p <- p +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 3) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
      labels = c("TRUE" = "\u663e\u8457", "FALSE" = "\u4e0d\u663e\u8457"),
      name   = "\u7edf\u8ba1\u663e\u8457\u6027"
    ) +
    ggplot2::scale_x_continuous(breaks = plot_data$n_pre) +
    ggplot2::labs(
      x        = "\u524d\u5904\u7406\u65f6\u671f\u6570\u91cf",
      y        = "ATT\u4f30\u8ba1",
      title    = "\u524d\u671f\u7a33\u5065\u6027\u89c4\u683c\u66f2\u7ebf",
      subtitle = sprintf("SR=%.1f%% (%s)", x$sensitivity_ratio * 100, x$robustness_level)
    ) +
    ggplot2::theme_minimal()

  p
}


# ---- Internal: anticipation sensitivity plot ----

.plot_anticipation_sensitivity <- function(x, ci_level = 0.95, show_ci = TRUE, ...) {
  # --- 1. Data preparation ---

  # Filter converged estimates
  converged_mask <- vapply(x$estimates, function(e) e$converged, logical(1))
  estimates <- x$estimates[converged_mask]

  # Empty data guard
  if (length(estimates) == 0L) {
    warning("\u6240\u6709\u4f30\u8ba1\u5747\u672a\u6536\u655b\uff0c\u65e0\u6cd5\u7ed8\u5236\u6709\u6548\u56fe\u5f62")
    p <- ggplot2::ggplot() +
      ggplot2::labs(
        title    = "\u65e0\u9884\u671f\u6548\u5e94\u654f\u611f\u6027\u5206\u6790",
        subtitle = "\u65e0\u6536\u655b\u4f30\u8ba1"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_text(color = "red")
      )
    return(p)
  }

  # Build plot_data data.frame using vapply
  plot_data <- data.frame(
    excluded = vapply(estimates, function(e) e$excluded_periods, integer(1)),
    att      = vapply(estimates, function(e) e$att, numeric(1)),
    ci_lower = vapply(estimates, function(e) e$ci_lower, numeric(1)),
    ci_upper = vapply(estimates, function(e) e$ci_upper, numeric(1))
  )

  baseline_att <- x$baseline_estimate$att

  # --- 2. Build ggplot ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = excluded, y = att))

  # Baseline ATT dashed line
  p <- p + ggplot2::geom_hline(yintercept = baseline_att, linetype = "dashed", color = "gray50")

  # Zero reference line
  p <- p + ggplot2::geom_hline(yintercept = 0, color = "gray80")

  # Confidence intervals
  if (show_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2
    )
  }

  # Point estimates and connection line
  p <- p +
    ggplot2::geom_point(size = 3, color = "#2166AC") +
    ggplot2::geom_line(color = "#2166AC", alpha = 0.5)

  # x-axis integer ticks
  p <- p + ggplot2::scale_x_continuous(breaks = plot_data$excluded)

  # --- 3. Recommended exclusion marker ---
  if (isTRUE(x$anticipation_detected)) {
    text_y <- if (show_ci) max(plot_data$ci_upper) else max(plot_data$att)

    p <- p +
      ggplot2::geom_vline(
        xintercept = x$recommended_exclusion,
        linetype   = "dotted",
        color      = "red",
        linewidth  = 1
      ) +
      ggplot2::annotate(
        "text",
        x     = x$recommended_exclusion + 0.2,
        y     = text_y,
        label = sprintf("\u63a8\u8350\u6392\u9664=%d", x$recommended_exclusion),
        color = "red",
        hjust = 0,
        size  = 3.5
      )
  }

  # --- 4. Labels ---
  subtitle <- if (isTRUE(x$anticipation_detected)) {
    "\u68c0\u6d4b\u5230\u9884\u671f\u6548\u5e94"
  } else {
    "\u672a\u68c0\u6d4b\u5230\u9884\u671f\u6548\u5e94"
  }

  p <- p +
    ggplot2::labs(
      x        = "\u6392\u9664\u7684\u5904\u7406\u524d\u65f6\u671f\u6570",
      y        = "ATT\u4f30\u8ba1",
      title    = "\u65e0\u9884\u671f\u6548\u5e94\u654f\u611f\u6027\u5206\u6790",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()

  p
}


# ---- S3 method: plot.lwdid_sensitivity ----

#' @export
plot.lwdid_sensitivity <- function(x, type = NULL, ci_level = 0.95,
                                    show_threshold = TRUE, show_ci = TRUE, ...) {
  stopifnot(inherits(x, "lwdid_sensitivity"))

  # Type inference
  if (is.null(type)) {
    type <- x$type
  }

  # Dispatch
  switch(type,
    "pre_period"      = .plot_specification_curve(x, ci_level = ci_level,
                                                   show_threshold = show_threshold,
                                                   show_ci = show_ci, ...),
    "no_anticipation" = .plot_anticipation_sensitivity(x, ci_level = ci_level,
                                                       show_ci = show_ci, ...),
    stop(sprintf("\u672a\u77e5\u7684\u654f\u611f\u6027\u5206\u6790\u7c7b\u578b: '%s'", type))
  )
}


#' @export
plot.lwdid_sensitivity_comprehensive <- function(x, show_ci = TRUE, ...) {
  stopifnot(inherits(x, "lwdid_sensitivity_comprehensive"))

  # Check ggplot2 and patchwork availability
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("\u9700\u8981\u5b89\u88c5 ggplot2 \u5305: install.packages('ggplot2')")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("\u9700\u8981\u5b89\u88c5 patchwork \u5305: install.packages('patchwork')")
  }

  # Generate sub-plots
  p1 <- .plot_specification_curve(x$pre_period, show_ci = show_ci, ...)
  p2 <- .plot_anticipation_sensitivity(x$no_anticipation, show_ci = show_ci, ...)

  # Vertical stack using patchwork
  p1 / p2
}


# ---- Internal: empty trend plot helper ----

.empty_trend_plot <- function(title, subtitle = "No plottable diagnostics") {
  ggplot2::ggplot() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_minimal()
}


# ---- S3 method: plot.lwdid_parallel_trends ----

#' @export
plot.lwdid_parallel_trends <- function(x, show_ci = TRUE, ...) {
  stopifnot(inherits(x, "lwdid_parallel_trends"))

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Install it with install.packages('ggplot2')")
  }

  estimates <- x$pre_trend_estimates %||% list()
  if (length(estimates) == 0L) {
    return(.empty_trend_plot(
      title = "Parallel Trends Diagnostics",
      subtitle = "No valid pre-treatment estimates"
    ))
  }

  plot_data <- do.call(rbind, lapply(estimates, function(estimate) {
    df_value <- max(as.integer(estimate$df %||% 1L), 1L)
    se_value <- as.numeric(estimate$se %||% NA_real_)
    att_value <- as.numeric(estimate$att %||% NA_real_)
    t_crit <- stats::qt(0.975, df = df_value)

    data.frame(
      event_time = as.integer(estimate$event_time %||% NA_integer_),
      att = att_value,
      se = se_value,
      ci_lower = att_value - t_crit * se_value,
      ci_upper = att_value + t_crit * se_value,
      pvalue = as.numeric(estimate$pvalue %||% NA_real_),
      stringsAsFactors = FALSE
    )
  }))

  subtitle <- if (length(x$joint_df %||% integer(0)) == 2L) {
    sprintf(
      "Joint F(%d, %d) = %s, p = %s",
      as.integer(x$joint_df[[1L]]),
      as.integer(x$joint_df[[2L]]),
      .format_trend_numeric(x$joint_f_stat),
      .format_pvalue(x$joint_pvalue, 4L)
    )
  } else {
    "Joint F-test unavailable"
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = event_time, y = att)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_line(color = "#2166AC", alpha = 0.6) +
    ggplot2::geom_point(
      ggplot2::aes(color = pvalue < 0.05),
      size = 2.5
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
      guide = "none"
    ) +
    ggplot2::labs(
      x = "Event time",
      y = "Placebo ATT",
      title = "Parallel Trends Diagnostics",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(show_ci)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.15,
      color = "gray40"
    )
  }

  p
}


# ---- Public function: plot_cohort_trends ----

#' Plot outcome trajectories by treatment cohort
#'
#' @param data data.frame or data.table in long format.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1). Cohort variable.
#' @param controls Optional character vector of controls (reserved for future use).
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @param normalize Logical scalar. Whether to subtract a baseline period.
#' @param normalize_period Optional baseline period.
#' @param show_treatment_lines Logical scalar.
#' @param show_trend_lines Logical scalar.
#' @param confidence_bands Logical scalar.
#' @param alpha Numeric significance level for confidence bands.
#' @return A ggplot object.
#' @export
plot_cohort_trends <- function(
    data,
    y,
    ivar,
    tvar,
    gvar,
    controls = NULL,
    never_treated_values = c(0, Inf),
    normalize = TRUE,
    normalize_period = NULL,
    show_treatment_lines = TRUE,
    show_trend_lines = TRUE,
    confidence_bands = TRUE,
    alpha = 0.05
) {
  required_columns <- c(y, ivar, tvar, gvar, controls %||% character(0))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Install it with install.packages('ggplot2')")
  }

  dt <- data.table::as.data.table(data.table::copy(data))
  custom_never <- setdiff(never_treated_values, c(0, Inf))
  if (length(custom_never) > 0L) {
    custom_mask <- dt[[gvar]] %in% custom_never
    custom_mask[is.na(custom_mask)] <- FALSE
    if (any(custom_mask)) {
      data.table::set(dt, i = which(custom_mask), j = gvar, value = Inf)
    }
  }

  cohorts <- .get_valid_cohorts(
    dt,
    gvar = gvar,
    ivar = ivar,
    never_treated_values = never_treated_values
  )
  if (length(cohorts) == 0L) {
    stop("No valid treatment cohorts found for plotting.")
  }

  t_min <- suppressWarnings(min(dt[[tvar]], na.rm = TRUE))
  t_max <- suppressWarnings(max(dt[[tvar]], na.rm = TRUE))
  if (isTRUE(normalize) && is.null(normalize_period)) {
    normalize_period <- min(cohorts) - 1L
  }

  cohort_plot_data <- lapply(cohorts, function(cohort_value) {
    cohort_mask <- abs(dt[[gvar]] - cohort_value) < LWDID_COHORT_FLOAT_TOLERANCE
    cohort_mask[is.na(cohort_mask)] <- FALSE
    cohort_data <- dt[cohort_mask]
    summary_data <- cohort_data[
      ,
      .(
        mean_y = mean(.SD[[1L]], na.rm = TRUE),
        sd_y = stats::sd(.SD[[1L]], na.rm = TRUE),
        n_y = sum(is.finite(.SD[[1L]]))
      ),
      by = tvar,
      .SDcols = y
    ]
    data.table::setorderv(summary_data, tvar)
    summary_data[, time := .SD[[1L]], .SDcols = tvar]

    if (isTRUE(normalize) && !is.null(normalize_period)) {
      baseline <- summary_data[time == normalize_period, mean_y]
      if (length(baseline) == 1L && is.finite(baseline)) {
        summary_data[, mean_y := mean_y - baseline]
      }
    }

    summary_data[, se_y := ifelse(n_y > 0, sd_y / sqrt(n_y), NA_real_)]
    summary_data[, t_crit := stats::qt(1 - alpha / 2, df = pmax(n_y - 1L, 1L))]
    summary_data[, ci_lower := mean_y - t_crit * se_y]
    summary_data[, ci_upper := mean_y + t_crit * se_y]
    summary_data[, cohort := as.factor(cohort_value)]
    summary_data[, treatment_time := cohort_value]
    summary_data[, is_control := FALSE]
    summary_data[]
  })

  plot_data <- data.table::rbindlist(cohort_plot_data, fill = TRUE)

  control_mask <- is_never_treated(dt[[gvar]])
  control_mask[is.na(control_mask)] <- FALSE
  control_data <- dt[control_mask]
  if (nrow(control_data) > 0L) {
    control_summary <- control_data[
      ,
      .(mean_y = mean(.SD[[1L]], na.rm = TRUE)),
      by = tvar,
      .SDcols = y
    ]
    data.table::setorderv(control_summary, tvar)
    control_summary[, time := .SD[[1L]], .SDcols = tvar]

    if (isTRUE(normalize) && !is.null(normalize_period)) {
      baseline <- control_summary[time == normalize_period, mean_y]
      if (length(baseline) == 1L && is.finite(baseline)) {
        control_summary[, mean_y := mean_y - baseline]
      }
    }

    control_summary[, cohort := as.factor("Never Treated")]
    control_summary[, treatment_time := NA_real_]
    control_summary[, is_control := TRUE]
    plot_data <- data.table::rbindlist(list(plot_data, control_summary), fill = TRUE)
  }

  p <- ggplot2::ggplot(
    plot_data[is_control == FALSE],
    ggplot2::aes(x = time, y = mean_y, color = cohort, fill = cohort)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "Time period",
      y = if (isTRUE(normalize)) "Outcome (normalized)" else "Outcome",
      title = "Outcome Trends by Treatment Cohort"
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(confidence_bands) && all(c("ci_lower", "ci_upper") %in% names(plot_data))) {
    p <- p + ggplot2::geom_ribbon(
      data = plot_data[is_control == FALSE],
      ggplot2::aes(
        x = time,
        ymin = ci_lower,
        ymax = ci_upper
      ),
      alpha = 0.15,
      color = NA,
      inherit.aes = FALSE
    )
  }

  if (isTRUE(show_trend_lines)) {
    trend_lines <- lapply(cohorts, function(cohort_value) {
      cohort_summary <- plot_data[
        is_control == FALSE & cohort == as.factor(cohort_value)
      ]
      pre_summary <- cohort_summary[time < cohort_value]
      if (nrow(pre_summary) < 2L) {
        return(NULL)
      }

      fit <- stats::lm(mean_y ~ time, data = pre_summary)
      trend_times <- seq.int(t_min, min(cohort_value + 1L, t_max))
      trend_values <- stats::predict(
        fit,
        newdata = data.frame(time = trend_times)
      )
      data.frame(
        time = trend_times,
        fitted = as.numeric(trend_values),
        cohort = as.factor(cohort_value),
        stringsAsFactors = FALSE
      )
    })
    trend_lines <- Filter(Negate(is.null), trend_lines)
    if (length(trend_lines) > 0L) {
      trend_df <- data.table::rbindlist(trend_lines, fill = TRUE)
      p <- p + ggplot2::geom_line(
        data = trend_df,
        ggplot2::aes(x = time, y = fitted, color = cohort),
        linetype = "dashed",
        alpha = 0.6,
        inherit.aes = FALSE
      )
    }
  }

  if (isTRUE(show_treatment_lines)) {
    for (cohort_value in cohorts) {
      p <- p + ggplot2::geom_vline(
        xintercept = cohort_value - 0.5,
        linetype = "dotted",
        color = "gray50",
        alpha = 0.5
      )
    }
  }

  if (nrow(control_data) > 0L) {
    p <- p + ggplot2::geom_line(
      data = plot_data[is_control == TRUE],
      ggplot2::aes(x = time, y = mean_y),
      color = "gray40",
      linetype = "dotdash",
      linewidth = 0.8,
      inherit.aes = FALSE
    )
  }

  p
}


# ---- S3 method: plot.lwdid_heterogeneous_trends ----

#' @export
plot.lwdid_heterogeneous_trends <- function(x, show_control = TRUE,
                                             show_ci = TRUE,
                                             facet_by_cohort = FALSE,
                                             show_slope_labels = TRUE, ...) {
  # --- 1. Input validation ---
  stopifnot(inherits(x, "lwdid_heterogeneous_trends"))

  # --- 2. Data preparation: build plot_data for each cohort trend ---
  plot_data <- data.frame()
  for (trend in x$trend_by_cohort) {
    cohort_g <- trend$cohort
    n_pre <- trend$n_pre_periods
    time_seq <- seq(cohort_g - n_pre, cohort_g - 1)
    t_mean <- mean(time_seq)
    t_centered <- time_seq - t_mean
    fitted_vals <- trend$intercept + trend$slope * t_centered
    t_crit <- stats::qt(0.975, df = max(n_pre - 2, 1))
    se_fitted <- sqrt(trend$slope_se^2 * t_centered^2)
    ci_lo <- fitted_vals - t_crit * se_fitted
    ci_hi <- fitted_vals + t_crit * se_fitted

    cohort_df <- data.frame(
      time = time_seq, fitted = fitted_vals,
      ci_lower = ci_lo, ci_upper = ci_hi,
      cohort = as.factor(cohort_g),
      treatment_time = cohort_g, is_control = FALSE
    )
    plot_data <- rbind(plot_data, cohort_df)
  }

  if (nrow(plot_data) == 0L) {
    return(.empty_trend_plot(
      title = "Cohort Trend Diagnostics",
      subtitle = "No cohort trends available"
    ))
  }

  # --- 3. Control group trend (optional) ---
  if (show_control && !is.null(x$control_group_trend)) {
    ctrl <- x$control_group_trend
    all_times <- sort(unique(plot_data$time))
    t_mean_ctrl <- mean(all_times)
    t_centered_ctrl <- all_times - t_mean_ctrl
    fitted_ctrl <- ctrl$intercept + ctrl$slope * t_centered_ctrl
    t_crit_ctrl <- stats::qt(0.975, df = max(length(all_times) - 2, 1))
    se_fitted_ctrl <- sqrt(ctrl$slope_se^2 * t_centered_ctrl^2)

    ctrl_df <- data.frame(
      time = all_times, fitted = fitted_ctrl,
      ci_lower = fitted_ctrl - t_crit_ctrl * se_fitted_ctrl,
      ci_upper = fitted_ctrl + t_crit_ctrl * se_fitted_ctrl,
      cohort = as.factor("\u5bf9\u7167\u7ec4"),
      treatment_time = NA, is_control = TRUE
    )
    plot_data <- rbind(plot_data, ctrl_df)
  }

  # --- 4. Build ggplot ---
  treated_data <- plot_data[!plot_data$is_control, ]
  p <- ggplot2::ggplot(treated_data,
                        ggplot2::aes(x = time, y = fitted,
                                     color = cohort, fill = cohort)) +
    ggplot2::geom_line(linewidth = 0.8)

  if (show_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      alpha = 0.15, color = NA
    )
  }

  # --- 5. Control group trend line (dashed) ---
  if (show_control && any(plot_data$is_control)) {
    ctrl_data <- plot_data[plot_data$is_control, ]
    p <- p + ggplot2::geom_line(
      data = ctrl_data,
      ggplot2::aes(x = time, y = fitted),
      linetype = "dashed", color = "gray30", linewidth = 0.8,
      inherit.aes = FALSE
    )
    if (show_ci) {
      p <- p + ggplot2::geom_ribbon(
        data = ctrl_data,
        ggplot2::aes(x = time, ymin = ci_lower, ymax = ci_upper),
        fill = "gray50", alpha = 0.1, inherit.aes = FALSE
      )
    }
  }

  # --- 6. Treatment time vertical lines ---
  treatment_times <- unique(stats::na.omit(plot_data$treatment_time))
  for (g_time in treatment_times) {
    p <- p + ggplot2::geom_vline(
      xintercept = g_time - 0.5,
      linetype = "dotted", color = "gray40", alpha = 0.7
    )
  }

  # --- 7. Slope labels (optional) ---
  if (show_slope_labels) {
    slope_labels <- data.frame()
    for (trend in x$trend_by_cohort) {
      g_val <- trend$cohort
      n_pre_val <- trend$n_pre_periods
      time_seq_sl <- seq(g_val - n_pre_val, g_val - 1)
      t_mean_sl <- mean(time_seq_sl)
      mid_time <- t_mean_sl
      mid_fitted <- trend$intercept
      p_str <- .format_pvalue(trend$slope_pvalue, 3L)
      slope_labels <- rbind(slope_labels, data.frame(
        time = mid_time, fitted = mid_fitted,
        cohort = as.factor(g_val),
        slope_label = sprintf("\u03b2=%.4f (p=%s)", trend$slope, p_str),
        slope_significant = (trend$slope_pvalue < 0.05),
        stringsAsFactors = FALSE
      ))
    }
    p <- p + ggplot2::geom_label(
      data = slope_labels,
      ggplot2::aes(x = time, y = fitted, label = slope_label,
                   color = cohort),
      fill = "white", alpha = 0.85, size = 3,
      fontface = ifelse(slope_labels$slope_significant, "bold", "plain"),
      label.padding = ggplot2::unit(0.15, "lines"),
      linewidth = 0.3,
      show.legend = FALSE, inherit.aes = FALSE
    )
  }

  # --- 8. Heterogeneity test subtitle ---
  het_label <- NULL
  if (!is.null(x$trend_heterogeneity_test)) {
    het_label <- sprintf(
      "\u5f02\u8d28\u6027F\u68c0\u9a8c: F(%d,%d)=%.2f, p=%s",
      x$trend_heterogeneity_test$df_num,
      x$trend_heterogeneity_test$df_den,
      x$trend_heterogeneity_test$f_stat,
      .format_pvalue(x$trend_heterogeneity_test$pvalue, 4L)
    )
  }

  # --- 9. Labels and theme ---
  p <- p +
    ggplot2::labs(
      x = "\u65f6\u671f", y = "\u62df\u5408\u7ed3\u679c\u53d8\u91cf",
      title = "\u961f\u5217\u5f02\u8d28\u6027\u8d8b\u52bf\u8bca\u65ad\u56fe",
      subtitle = het_label,
      color = "\u961f\u5217", fill = "\u961f\u5217"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  # --- 10. Faceting ---
  if (facet_by_cohort) {
    p <- p + ggplot2::facet_wrap(~cohort, scales = "free_y")
  }

  # --- 11. Return ---
  p
}

# ---- S3 method: plot.lwdid_transformation_recommendation ----

#' @export
plot.lwdid_transformation_recommendation <- function(x, ...) {
  # --- 1. Input validation ---
  stopifnot(inherits(x, "lwdid_transformation_recommendation"))

  # --- 2. Data preparation ---
  df <- data.frame(
    method = names(x$scores),
    score = as.numeric(x$scores),
    is_recommended = (names(x$scores) == x$recommended_method),
    stringsAsFactors = FALSE
  )
  # Order by score descending, rev so highest is at top in coord_flip
  df$method <- factor(df$method,
                       levels = rev(df$method[order(df$score)]))

  # Method labels replacement (partial coverage OK)
  if (!is.null(x$method_labels)) {
    new_levels <- x$method_labels[levels(df$method)]
    has_label <- !is.na(new_levels)
    final_levels <- levels(df$method)
    final_levels[has_label] <- new_levels[has_label]
    levels(df$method) <- final_levels
  }

  # --- 3. Determine plot mode ---

  if (!is.null(x$sub_scores)) {
    # Mode A: stacked bar chart (sub_scores available)
    sub_df <- do.call(rbind, lapply(names(x$sub_scores), function(m) {
      ss <- x$sub_scores[[m]]
      data.frame(method = m, dimension = names(ss),
                 sub_score = as.numeric(ss),
                 stringsAsFactors = FALSE)
    }))
    sub_df$method <- factor(sub_df$method,
                             levels = levels(df$method))
    sub_df$is_recommended <- (sub_df$method ==
                                x$recommended_method)
    sub_df$dimension <- factor(sub_df$dimension,
                                levels = unique(sub_df$dimension))

    p <- ggplot2::ggplot(
      sub_df,
      ggplot2::aes(x = method, y = sub_score,
                    fill = dimension)
    ) +
      ggplot2::geom_col(
        ggplot2::aes(alpha = is_recommended),
        width = 0.6, position = "stack"
      ) +
      ggplot2::scale_alpha_manual(
        values = c("TRUE" = 1.0, "FALSE" = 0.6),
        guide = "none"
      ) +
      # Recommended method border highlight
      ggplot2::geom_col(
        data = sub_df[sub_df$is_recommended, ],
        ggplot2::aes(x = method, y = sub_score),
        fill = NA, color = "#2166AC", linewidth = 1.0,
        width = 0.6, position = "stack",
        show.legend = FALSE
      ) +
      # Total score annotation at bar top
      ggplot2::geom_text(
        data = df,
        ggplot2::aes(x = method, y = score,
                      label = sprintf("%.2f", score)),
        hjust = -0.3, size = 3.5, inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_brewer(
        palette = "Set2",
        name = "\u8bc4\u5206\u7ef4\u5ea6"
      )
  } else {
    # Mode B: simple bar chart (no sub_scores)
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = method, y = score,
                    fill = is_recommended)
    ) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#2166AC", "FALSE" = "gray60"),
        guide = "none"
      ) +
      # Recommended method border highlight
      ggplot2::geom_col(
        data = df[df$is_recommended, ],
        ggplot2::aes(x = method, y = score),
        fill = NA, color = "#2166AC", linewidth = 1.2,
        width = 0.6, show.legend = FALSE
      ) +
      # Score annotation
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", score)),
        hjust = -0.3, size = 3.5
      )
  }

  # --- 4. Labels and theme ---

  # Confidence level color mapping
  conf_color <- switch(x$confidence_level,
    "High"   = "#2166AC",
    "Medium" = "#F4A582",
    "Low"    = "#B2182B",
    "gray50"
  )

  p <- p +
    ggplot2::labs(
      x = NULL,
      y = "\u7efc\u5408\u8bc4\u5206",
      title = "\u53d8\u6362\u65b9\u6cd5\u63a8\u8350\u8bc4\u5206",
      subtitle = sprintf(
        "\u63a8\u8350\u65b9\u6cd5: %s\uff08\u7f6e\u4fe1\u5ea6: %s\uff09",
        x$recommended_method, x$confidence_level
      )
    ) +
    ggplot2::coord_flip(
      ylim = c(0, max(df$score) * 1.2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (!is.null(x$sub_scores))
        "bottom" else "none",
      plot.subtitle = ggplot2::element_text(
        color = conf_color
      )
    )

  p
}


# ============================================================================
# Internal diagnostic panel helpers (used by plot.R .plot_diagnostics_panel)
# ============================================================================

#' @keywords internal
.plot_trends_diagnostic <- function(x) {
  # x is a lwdid_trends object with pre_treatment_coefficients and f_test
  df <- x$pre_treatment_coefficients
  if (is.null(df) || nrow(df) == 0L) {
    return(ggplot2::ggplot() +
      ggplot2::labs(title = "Parallel Trends: No data") +
      ggplot2::theme_minimal())
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = period, y = coefficient)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "gray50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2, color = "gray40"
    ) +
    ggplot2::geom_point(size = 3, color = "#2166AC") +
    ggplot2::geom_line(color = "#2166AC", alpha = 0.5)

  # F-test subtitle
  f_test <- x$f_test
  subtitle <- NULL
  if (!is.null(f_test)) {
    subtitle <- sprintf("F(%d,%d)=%.2f, p=%s",
      f_test$df1, f_test$df2, f_test$f_stat,
      .format_pvalue(f_test$f_pvalue))
  }

  p <- p +
    ggplot2::labs(
      x = "Pre-treatment Period",
      y = "Coefficient",
      title = "Parallel Trends Test",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  p
}

#' @keywords internal
.plot_clustering_diagnostic <- function(x) {
  # x is a lwdid_clustering_diagnosis object
  sizes <- x$cluster_sizes
  if (is.null(sizes) || length(sizes) == 0L) {
    return(ggplot2::ggplot() +
      ggplot2::labs(title = "Clustering: No data") +
      ggplot2::theme_minimal())
  }

  df <- data.frame(
    cluster = names(sizes),
    size = as.integer(sizes),
    stringsAsFactors = FALSE
  )
  # Sort by size descending
  df <- df[order(-df$size), ]

  # Truncate for large cluster counts (>60 → top/bottom 15)
  n_clusters <- nrow(df)
  suffix <- ""
  if (n_clusters > 60L) {
    top <- utils::head(df, 15L)
    bottom <- utils::tail(df, 15L)
    df <- rbind(top, bottom)
    suffix <- sprintf(" (showing top/bottom 15 of %d)", n_clusters)
  }

  df$cluster <- factor(df$cluster, levels = df$cluster)

  # Mean reference line value
  mean_size <- mean(as.integer(sizes))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = size)) +
    ggplot2::geom_col(fill = "#2166AC", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = mean_size, linetype = "dashed",
                        color = "red", alpha = 0.7) +
    ggplot2::labs(
      x = "Cluster",
      y = "Size",
      title = "Cluster Size Distribution",
      subtitle = sprintf("N=%d, Balance=%.1f, ICC=%.3f, Eff.clusters=%.1f%s",
        x$n_clusters, x$balance_ratio, x$icc,
        x$effective_clusters, suffix)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7)
    )

  # >30 clusters: use horizontal bars (coord_flip) to avoid label overlap
  if (n_clusters > 30L) {
    p <- p + ggplot2::coord_flip()
  }

  p
}

#' @keywords internal
.plot_selection_diagnostic <- function(x) {
  # x is a lwdid_selection_diagnosis object
  df <- .coerce_selection_attrition_series(x$attrition_by_period)
  if (is.null(df) || nrow(df) == 0L) {
    # Text summary fallback
    summary_text <- sprintf(
      "Overall Attrition: %.1f%%\nRisk Level: %s",
      (x$attrition_rate %||% 0) * 100,
      x$selection_risk %||% "Unknown"
    )
    p <- ggplot2::ggplot(data.frame(x = 0.5, y = 0.5, label = summary_text),
                          ggplot2::aes(x = x, y = y, label = label)) +
      ggplot2::geom_text(size = 5) +
      ggplot2::labs(title = "Selection / Attrition Summary") +
      ggplot2::theme_void()
    return(p)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = period, y = attrition_rate)) +
    ggplot2::geom_line(color = "#B2182B", linewidth = 1) +
    ggplot2::geom_point(color = "#B2182B", size = 2) +
    ggplot2::labs(
      x = "Period",
      y = "Attrition Rate",
      title = "Selection / Attrition Over Time",
      subtitle = sprintf("Overall attrition: %.1f%%, Risk: %s",
        x$attrition_rate * 100, x$selection_risk)
    ) +
    ggplot2::theme_minimal()

  # Format y-axis as percentage
  has_scales <- requireNamespace("scales", quietly = TRUE)
  if (has_scales) {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1))
  }
  p
}

#' @keywords internal
.coerce_selection_attrition_series <- function(attrition_by_period) {
  if (is.null(attrition_by_period)) {
    return(NULL)
  }

  if (is.data.frame(attrition_by_period)) {
    return(attrition_by_period)
  }

  values <- unlist(attrition_by_period, use.names = FALSE)
  if (length(values) == 0L) {
    return(NULL)
  }

  periods <- names(attrition_by_period)
  if (is.null(periods)) {
    periods <- seq_along(values)
  } else {
    periods <- utils::type.convert(periods, as.is = TRUE)
  }

  data.frame(
    period = periods,
    attrition_rate = as.numeric(values),
    stringsAsFactors = FALSE
  )
}
