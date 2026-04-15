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
    stop("\u9700\u8981\u5b89\u88c5ggplot2\u5305\u3002\u8fd0\u884c: install.packages('ggplot2')",
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
        stop("\u65e0\u654f\u611f\u6027\u5206\u6790\u7ed3\u679c\u3002\u8bf7\u5148\u8fd0\u884csensitivity_analysis()",
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
    stop("\u65e0\u7ed8\u56fe\u6570\u636e\u3002\u8bf7\u5728lwdid()\u4e2d\u8bbe\u7f6egraph=TRUE",
         call. = FALSE)

  plot_data <- x$metadata$plot_data

  # --- Staggered mode: gid dispatch ---
  if (isTRUE(x$is_staggered)) {
    if (is.null(plot_data$cohorts))
      stop("Staggered\u6a21\u5f0f\u4e0bplot_data\u7f3a\u5c11cohorts\u5b57\u6bb5",
           call. = FALSE)

    if (!is.null(gid)) {
      # Specific cohort
      gid_str <- as.character(gid)
      if (!(gid_str %in% names(plot_data$cohorts)))
        stop(sprintf(
          "gid=%s \u4e0d\u5728\u53ef\u7528cohort\u4e2d\u3002\u53ef\u7528: %s",
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
      warning("gid\u53c2\u6570\u4ec5\u5728Staggered\u6a21\u5f0f\u4e0b\u751f\u6548\uff0c\u5df2\u5ffd\u7565",
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
    stop("\u65e0\u8bca\u65ad\u4fe1\u606f\u3002\u8bf7\u5728lwdid()\u4e2d\u8bbe\u7f6ereturn_diagnostics=TRUE",
         call. = FALSE)

  # Validate which parameter
  valid_which <- c("parallel_trends", "clustering", "selection")
  if (!is.null(which)) {
    which <- match.arg(which, valid_which, several.ok = TRUE)
  }

  should_include <- function(name) is.null(which) || name %in% which

  panels <- list()

  if (should_include("parallel_trends") &&
      !is.null(x$diagnostics$parallel_trends)) {
    panels$trends <- .plot_trends_diagnostic(
      x$diagnostics$parallel_trends)
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
      available <- intersect(
        names(x$diagnostics),
        valid_which[vapply(valid_which, function(w)
          !is.null(x$diagnostics[[w]]), logical(1))]
      )
      stop(sprintf(
        "\u6307\u5b9a\u7684\u8bca\u65ad\u7c7b\u578b [%s] \u65e0\u53ef\u7528\u6570\u636e\u3002\u53ef\u7528\u8bca\u65ad: %s",
        paste(which, collapse = ", "),
        if (length(available) > 0)
          paste(available, collapse = ", ") else "\u65e0"
      ), call. = FALSE)
    } else {
      stop("\u8bca\u65ad\u5bf9\u8c61\u4e2d\u65e0\u53ef\u7ed8\u5236\u5185\u5bb9\uff08\u6240\u6709\u8bca\u65ad\u5b50\u9879\u5747\u4e3aNULL\uff09",
           call. = FALSE)
    }
  }

  if (length(panels) == 1L) return(panels[[1L]])

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("\u9700\u8981\u5b89\u88c5patchwork\u5305\u7528\u4e8e\u9762\u677f\u56fe\u3002\u8fd0\u884c: install.packages('patchwork')",
         call. = FALSE)

  combined <- Reduce(`+`, panels) +
    patchwork::plot_layout(ncol = min(2L, length(panels))) +
    patchwork::plot_annotation(title = "Diagnostics Panel")
  combined
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
    stop("\u9700\u8981\u5b89\u88c5ggplot2\u5305", call. = FALSE)
  stopifnot(inherits(x, "lwdid_trends"))
  type <- match.arg(type, c("coefficients", "trajectories"))

  switch(type,
    "coefficients" = .plot_trends_diagnostic(x),
    "trajectories" = {
      if (is.null(x$group_trends))
        stop("\u65e0\u5206\u7ec4\u8d8b\u52bf\u6570\u636e", call. = FALSE)

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
            "F(%d,%d)=%.2f, p=%s | smooth=%s",
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
    stop("\u9700\u8981\u5b89\u88c5ggplot2\u5305", call. = FALSE)
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
    stop("\u9700\u8981\u5b89\u88c5ggplot2\u5305", call. = FALSE)
  stopifnot(inherits(x, "lwdid_selection_diagnosis"))
  .plot_selection_diagnostic(x)
}
