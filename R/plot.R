# ===========================================================================
# plot.R — Unified plotting entry point and diagnostic object plot methods
#
# Implements plot.lwdid_result() as unified entry, dispatching to different
# visualization functions by type. Supports four plot types:
# event_study, trajectories, sensitivity, diagnostics.
# Also implements S3 plot methods for diagnostic objects.
# ===========================================================================

#' @title Plot lwdid results
#' @description Unified plotting entry point. Dispatches to specific
#'   visualization functions based on the \code{type} parameter.
#' @param x lwdid_result object
#' @param type character, plot type: "event_study" (default),
#'   "trajectories", "sensitivity", "diagnostics"
#' @param ... passed to the specific plotting function
#' @return ggplot object
#' @export
plot.lwdid_result <- function(x, type = "event_study", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting. Run install.packages('ggplot2').",
         call. = FALSE)
  }
  stopifnot(inherits(x, "lwdid_result"))
  type <- match.arg(type, c("event_study", "trajectories",
                              "sensitivity", "diagnostics"))

  switch(type,
    "event_study"  = plot_event_study(x, ...),
    "trajectories" = .plot_trajectories(x, ...),
    "sensitivity"  = {
      if (is.null(x$sensitivity))
        stop("No sensitivity-analysis results are available. Run sensitivity_analysis() first.",
             call. = FALSE)
      plot(x$sensitivity, ...)
    },
    "diagnostics"  = .plot_diagnostics_panel(x, ...)
  )
}


# ── .plot_trajectories ──────────────────────────────────────────────────────

#' @keywords internal
.plot_trajectories <- function(x, gid = NULL, show_individual = FALSE,
                                highlight_id = NULL, title = NULL,
                                theme = ggplot2::theme_minimal, ...) {
  if (is.null(x$metadata$plot_data))
    stop("No plotting data are available. Set graph = TRUE when calling lwdid().",
         call. = FALSE)

  plot_data <- x$metadata$plot_data

  # --- Staggered mode: gid dispatch ---
  if (isTRUE(x$is_staggered)) {
    if (is.null(plot_data$cohorts))
      stop("Staggered plotting data must include a cohorts field.",
           call. = FALSE)

    if (!is.null(gid)) {
      # Specific cohort
      gid_str <- as.character(gid)
      if (!(gid_str %in% names(plot_data$cohorts)))
        stop(sprintf(
          "gid=%s is not an available cohort. Available cohorts: %s",
          gid_str,
          paste(names(plot_data$cohorts), collapse = ", ")),
          call. = FALSE)
      pd <- plot_data$cohorts[[gid_str]]
      effective_title <- sprintf("Cohort %s: Treated vs Control",
                                 gid_str)
    } else {
      # gid=NULL: faceted plot of all cohorts
      all_dfs <- lapply(names(plot_data$cohorts), function(g) {
        cdata <- plot_data$cohorts[[g]]
        rbind(
          data.frame(time = cdata$time, y = cdata$control_mean,
                     group = "Control", cohort = g),
          data.frame(time = cdata$time, y = cdata$treated_mean,
                     group = "Treated", cohort = g)
        )
      })
      df_all <- do.call(rbind, all_dfs)

      # Per-cohort intervention lines
      intervention_df <- data.frame(
        cohort = names(plot_data$cohorts),
        xint = vapply(plot_data$cohorts,
                      function(c) c$intervention_point, numeric(1))
      )

      p <- ggplot2::ggplot(df_all, ggplot2::aes(
        x = time, y = y, color = group, linetype = group)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_vline(
          data = intervention_df,
          ggplot2::aes(xintercept = xint),
          linetype = "dashed", color = "black", alpha = 0.7,
          inherit.aes = FALSE
        ) +
        ggplot2::facet_wrap(~ cohort, scales = "free_x") +
        ggplot2::scale_color_manual(
          values = c("Control" = "#2166AC", "Treated" = "#B2182B")
        ) +
        ggplot2::scale_linetype_manual(
          values = c("Control" = "dashed", "Treated" = "solid")
        ) +
        ggplot2::labs(
          x = "Time", y = "Residualized Outcome",
          title = if (!is.null(title)) title
                  else "All Cohorts: Treated vs Control Trajectories",
          color = "Group", linetype = "Group"
        ) +
        theme()
      return(p)
    }
  } else {
    # Non-staggered: ignore gid
    if (!is.null(gid))
      warning("The gid argument only applies in staggered mode and has been ignored.",
              call. = FALSE)
    pd <- plot_data
    effective_title <- NULL
  }

  # --- Single cohort / non-staggered base plot ---
  df_ctrl <- data.frame(time = pd$time, y = pd$control_mean,
                         group = "Control")
  df_trt  <- data.frame(time = pd$time, y = pd$treated_mean,
                         group = "Treated")
  df <- rbind(df_ctrl, df_trt)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = y,
                                          color = group,
                                          linetype = group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = pd$intervention_point,
                         linetype = "dashed", color = "black",
                         alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c("Control" = "#2166AC", "Treated" = "#B2182B")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Control" = "dashed", "Treated" = "solid")
    ) +
    ggplot2::labs(
      x = "Time", y = "Residualized Outcome",
      title = if (!is.null(effective_title)) effective_title
              else if (!is.null(title)) title
              else "Treated vs Control Trajectories",
      color = "Group", linetype = "Group"
    ) +
    theme()

  # Individual trajectories overlay
  if (show_individual && !is.null(pd$treated_individual)) {
    df_indiv <- do.call(rbind, lapply(
      names(pd$treated_individual), function(id) {
        data.frame(time = pd$time,
                   y = pd$treated_individual[[id]],
                   unit_id = id)
      }))
    p <- p + ggplot2::geom_line(
      data = df_indiv,
      ggplot2::aes(x = time, y = y, group = unit_id),
      color = "gray70", alpha = 0.3, linewidth = 0.5,
      inherit.aes = FALSE
    )
  }

  # Highlight specific units
  if (!is.null(highlight_id) && !is.null(pd$treated_individual)) {
    valid_ids <- intersect(highlight_id,
                           names(pd$treated_individual))
    if (length(valid_ids) > 0L) {
      df_hl <- do.call(rbind, lapply(valid_ids, function(id) {
        data.frame(time = pd$time,
                   y = pd$treated_individual[[id]],
                   unit_id = id)
      }))
      p <- p + ggplot2::geom_line(
        data = df_hl,
        ggplot2::aes(x = time, y = y, group = unit_id),
        color = "#E66100", linewidth = 1.2,
        inherit.aes = FALSE
      )
    }
  }

  p
}


# ── .plot_diagnostics_panel ─────────────────────────────────────────────────

#' @keywords internal
.plot_diagnostics_panel <- function(x, which = NULL, ...) {
  if (is.null(x$diagnostics))
    stop("No diagnostics are available. Set return_diagnostics = TRUE when calling lwdid().",
         call. = FALSE)

  # Validate which parameter
  valid_which <- c("parallel_trends", "trends", "clustering", "selection")
  if (!is.null(which)) {
    which <- match.arg(which, valid_which, several.ok = TRUE)
    which[which == "trends"] <- "parallel_trends"
  }

  should_include <- function(name) is.null(which) || name %in% which

  panels <- list()

  trend_diagnostic <- x$diagnostics$parallel_trends %||% x$diagnostics$trends

  if (should_include("parallel_trends") &&
      !is.null(trend_diagnostic)) {
    panels$trends <- .plot_integrated_trend_diagnostic(trend_diagnostic)
  }
  if (should_include("clustering") &&
      !is.null(x$diagnostics$clustering)) {
    panels$clustering <- .plot_clustering_diagnostic(
      x$diagnostics$clustering)
  }
  if (should_include("selection") &&
      !is.null(x$diagnostics$selection)) {
    panels$selection <- .plot_selection_diagnostic(
      x$diagnostics$selection)
  }

  if (length(panels) == 0L) {
    if (!is.null(which)) {
      available_sources <- list(
        parallel_trends = trend_diagnostic,
        trends = trend_diagnostic,
        clustering = x$diagnostics$clustering,
        selection = x$diagnostics$selection
      )
      available <- names(available_sources)[vapply(
        available_sources, Negate(is.null), logical(1))]
      stop(sprintf(
        "The requested diagnostics [%s] have no available data. Available diagnostics: %s",
        paste(which, collapse = ", "),
        if (length(available) > 0)
          paste(available, collapse = ", ") else "none"
      ), call. = FALSE)
    } else {
      stop("The diagnostics object has no plottable content; all diagnostic entries are NULL.",
           call. = FALSE)
    }
  }

  if (length(panels) == 1L) return(panels[[1L]])

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("The patchwork package is required for diagnostics panel plots. Run install.packages('patchwork').",
         call. = FALSE)

  combined <- Reduce(`+`, panels) +
    patchwork::plot_layout(ncol = min(2L, length(panels))) +
    patchwork::plot_annotation(title = "Diagnostics Panel")
  combined
}

#' @keywords internal
.plot_integrated_trend_diagnostic <- function(x) {
  if (inherits(x, "lwdid_transformation_recommendation")) {
    return(plot.lwdid_transformation_recommendation(x))
  }

  if (inherits(x, "lwdid_parallel_trends")) {
    return(plot.lwdid_parallel_trends(x))
  }

  .plot_trends_diagnostic(x)
}


# ── plot.lwdid_trends ───────────────────────────────────────────────────────

#' @title Plot parallel trends diagnostic
#' @description Plot pre-treatment coefficient test or group trend
#'   trajectories for a \code{lwdid_trends} object.
#' @param x lwdid_trends object
#' @param type character, "coefficients" (default) or "trajectories"
#' @param smooth_method character, smoothing method for trajectories
#'   ("lm" default, or "loess")
#' @param ... passed to ggplot2
#' @return ggplot object
#' @export
plot.lwdid_trends <- function(x, type = "coefficients",
                               smooth_method = "lm", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  stopifnot(inherits(x, "lwdid_trends"))
  type <- match.arg(type, c("coefficients", "trajectories"))

  switch(type,
    "coefficients" = .plot_trends_diagnostic(x),
    "trajectories" = {
      if (is.null(x$group_trends))
        stop("No group-trend data are available.", call. = FALSE)

      df <- x$group_trends  # period, group, mean_y
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = period, y = mean_y,
        color = group, linetype = group)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_manual(
          values = c("treated" = "#B2182B",
                     "control" = "#2166AC")
        ) +
        ggplot2::scale_linetype_manual(
          values = c("treated" = "solid",
                     "control" = "dashed")
        ) +
        ggplot2::geom_smooth(
          method = smooth_method, se = TRUE,
          alpha = 0.15, linewidth = 0.5,
          formula = y ~ x) +
        ggplot2::labs(
          title = "Group Trend Trajectories",
          subtitle = sprintf(
            "Joint F-test: F(%d, %d) = %.2f, p = %s; smoother: %s",
            x$f_test$df1, x$f_test$df2,
            x$f_test$f_stat,
            .format_pvalue(x$f_test$f_pvalue),
            smooth_method),
          x = "Period", y = "Mean Outcome",
          color = "Group", linetype = "Group"
        ) +
        ggplot2::theme_minimal()
      p
    }
  )
}


# ── plot.lwdid_clustering_diagnosis ─────────────────────────────────────────

#' @title Plot clustering diagnostic
#' @description Plot cluster size distribution for a
#'   \code{lwdid_clustering_diagnosis} object.
#' @param x lwdid_clustering_diagnosis object
#' @param ... passed to ggplot2
#' @return ggplot object
#' @export
plot.lwdid_clustering_diagnosis <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  stopifnot(inherits(x, "lwdid_clustering_diagnosis"))
  .plot_clustering_diagnostic(x)
}


# ── plot.lwdid_selection_diagnosis ──────────────────────────────────────────

#' @title Plot selection diagnostic
#' @description Plot attrition over time for a
#'   \code{lwdid_selection_diagnosis} object.
#' @param x lwdid_selection_diagnosis object
#' @param ... passed to ggplot2
#' @return ggplot object
#' @export
plot.lwdid_selection_diagnosis <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  stopifnot(inherits(x, "lwdid_selection_diagnosis"))
  .plot_selection_diagnostic(x)
}
