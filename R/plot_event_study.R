# ============================================================================
# plot_event_study.R — Event Study visualization for lwdid
# ============================================================================
#
# Internal aggregation engine for event-study plots.
# Provides cohort-to-event-time aggregation with support for
# simple mean and cohort-size-weighted aggregation methods.
# ============================================================================


# ----------------------------------------------------------------------------
# .select_df_for_event_time
# ----------------------------------------------------------------------------
#' Select degrees of freedom for a given event time
#'
#' @param df_values Numeric vector of df values from individual cohort estimates.
#' @param weights   Optional numeric vector of weights (used when strategy = "weighted").
#' @param strategy  Character: "conservative", "weighted", or "fallback".
#' @param n_cohorts Integer number of cohorts contributing to this event time.
#' @return Integer degrees of freedom (>= 1).
#' @keywords internal
.select_df_for_event_time <- function(df_values, weights = NULL, strategy, n_cohorts) {
  # Filter valid df values: exclude NA and non-positive (<= 0)
  valid_mask <- !is.na(df_values) & df_values > 0
  valid_df_values <- df_values[valid_mask]

  # Fallback when all df values are invalid
  if (length(valid_df_values) == 0L) {
    return(max(1L, as.integer(n_cohorts - 1L)))
  }

  switch(strategy,
    "conservative" = {
      max(1L, as.integer(min(valid_df_values)))
    },
    "weighted" = {
      if (!is.null(weights)) {
        valid_weights <- weights[valid_mask]
        w_norm <- valid_weights / sum(valid_weights)
      } else {
        # aggregation = "mean" case: equal weights -> df(e) = round(mean(valid_df))
        w_norm <- rep(1 / length(valid_df_values), length(valid_df_values))
      }
      max(1L, as.integer(round(sum(w_norm * valid_df_values))))
    },
    "fallback" = {
      max(1L, as.integer(n_cohorts - 1L))
    },
    stop(sprintf("Unknown df_strategy: '%s'. Must be one of 'conservative', 'weighted', 'fallback'.", strategy))
  )
}


# ----------------------------------------------------------------------------
# .aggregate_by_event_time
# ----------------------------------------------------------------------------
#' Aggregate cohort-level ATT estimates to event-time level
#'
#' @param cohort_df   Data frame with columns: cohort, period, att, se, df_inference.
#'                    Column \code{event_time} is auto-calculated as \code{period - cohort}
#'                    if not already present.
#' @param method      Character: \code{"mean"} for simple average or \code{"weighted"}
#'                    for cohort-size-weighted aggregation.
#' @param cohort_sizes Named numeric vector mapping cohort identifiers (as character)
#'                     to sample sizes. Required when \code{method = "weighted"}.
#' @param df_strategy Character: strategy for aggregating degrees of freedom.
#'                    One of \code{"conservative"}, \code{"weighted"}, or \code{"fallback"}.
#' @return A \code{data.frame} with columns: event_time, att, se, df_inference, n_cohorts,
#'         sorted by event_time ascending.
#' @keywords internal
.aggregate_by_event_time <- function(cohort_df,
                                     method,
                                     cohort_sizes = NULL,
                                     df_strategy = "conservative") {

  # --------------------------------------------------------------------------
  # 1. event_time auto-calculation (if column doesn't exist)
  # --------------------------------------------------------------------------
  if (!"event_time" %in% names(cohort_df)) {
    required_cols <- c("cohort", "period")
    missing_cols <- setdiff(required_cols, names(cohort_df))
    if (length(missing_cols) > 0L) {
      stop(sprintf(
        paste0("att_by_cohort_time missing event_time column, ",
               "and cannot compute automatically from period-cohort ",
               "(missing columns: %s)"),
        paste(missing_cols, collapse = ", ")
      ))
    }
    cohort_df$event_time <- cohort_df$period - cohort_df$cohort
  }


  # --------------------------------------------------------------------------
  # 2. Iterate over unique event times and aggregate
  # --------------------------------------------------------------------------
  unique_event_times <- sort(unique(cohort_df$event_time))

  # Pre-allocate result lists
  res_event_time  <- numeric(0L)
  res_att         <- numeric(0L)
  res_se          <- numeric(0L)
  res_df          <- integer(0L)
  res_n_cohorts   <- integer(0L)

  for (e in unique_event_times) {
    et_data <- cohort_df[cohort_df$event_time == e, , drop = FALSE]

    # Remove rows where att or se is NA
    complete_mask <- !is.na(et_data$att) & !is.na(et_data$se)
    et_data <- et_data[complete_mask, , drop = FALSE]

    # Skip all-NA event_times
    if (nrow(et_data) == 0L) next

    G_e <- unique(et_data$cohort)
    n_g <- length(G_e)

    if (method == "mean") {
      att_e <- mean(et_data$att)
      se_e  <- sqrt(sum(et_data$se^2)) / n_g

      # df selection — no weights for mean method
      df_e <- .select_df_for_event_time(
        df_values = et_data$df_inference,
        weights   = NULL,
        strategy  = df_strategy,
        n_cohorts = n_g
      )

    } else if (method == "weighted") {
      # Look up cohort sizes; missing cohorts default to 0
      sizes <- cohort_sizes[as.character(G_e)]
      sizes[is.na(sizes)] <- 0

      total_size <- sum(sizes, na.rm = TRUE)
      if (total_size <= 0) {
        stop("All contributing cohorts for this event_time have size 0 or missing in cohort_sizes")
      }

      # Per-row weights aligned with et_data rows via cohort matching
      w_vec <- as.numeric(cohort_sizes[as.character(et_data$cohort)])
      w_vec[is.na(w_vec)] <- 0
      w_vec <- w_vec / total_size

      att_e <- sum(w_vec * et_data$att)
      se_e  <- sqrt(sum(w_vec^2 * et_data$se^2))

      # df selection — pass per-row weights
      df_e <- .select_df_for_event_time(
        df_values = et_data$df_inference,
        weights   = w_vec,
        strategy  = df_strategy,
        n_cohorts = n_g
      )

    } else {
      stop(sprintf("Unknown method: '%s'. Must be 'mean' or 'weighted'.", method))
    }

    # Accumulate results
    res_event_time <- c(res_event_time, e)
    res_att        <- c(res_att, att_e)
    res_se         <- c(res_se, se_e)
    res_df         <- c(res_df, df_e)
    res_n_cohorts  <- c(res_n_cohorts, as.integer(n_g))
  }


  # --------------------------------------------------------------------------
  # 3. Build output data.frame, sorted by event_time ascending
  # --------------------------------------------------------------------------
  result <- data.frame(
    event_time   = res_event_time,
    att          = res_att,
    se           = res_se,
    df_inference = res_df,
    n_cohorts    = res_n_cohorts,
    stringsAsFactors = FALSE
  )

  # Sort by event_time (should already be sorted, but enforce)
  result <- result[order(result$event_time), , drop = FALSE]
  rownames(result) <- NULL

  result
}


# ----------------------------------------------------------------------------
# .format_subtitle
# ----------------------------------------------------------------------------
#' Format subtitle string for event study plot
#'
#' @param aggregation Character: "mean" or "weighted".
#' @param df_strategy Character: "conservative", "weighted", or "fallback".
#' @param df_values   Numeric vector of df_inference values from plot data.
#' @param is_staggered Logical: whether the result is from staggered mode.
#' @return Character string for use as plot subtitle.
#' @keywords internal
.format_subtitle <- function(aggregation, df_strategy, df_values, is_staggered) {
  # Filter valid df values (exclude NA and Inf from anchors)
  valid_dfs <- df_values[is.finite(df_values) & !is.na(df_values)]

  if (!is_staggered) {
    if (length(valid_dfs) > 0L) {
      return(sprintf("Common Timing | df=%d", as.integer(valid_dfs[1L])))
    } else {
      return("Common Timing")
    }
  }

  # Staggered mode
  if (length(valid_dfs) == 0L) {
    return(sprintf("aggregation=%s | df_strategy=%s", aggregation, df_strategy))
  }

  df_min <- as.integer(min(valid_dfs))
  df_max <- as.integer(max(valid_dfs))

  if (df_min == df_max) {
    df_text <- sprintf("df=%d", df_min)
  } else {
    df_text <- sprintf("df range=[%d, %d]", df_min, df_max)
  }

  sprintf("aggregation=%s | df_strategy=%s | %s", aggregation, df_strategy, df_text)
}


# ----------------------------------------------------------------------------
# plot_event_study.lwdid_result — S3 method
# ----------------------------------------------------------------------------
#' Event Study Plot for lwdid Results
#'
#' @param x An \code{lwdid_result} object.
#' @param ci_level Numeric in (0,1): confidence interval level (default 0.95).
#' @param reference_line Numeric: y-value for horizontal reference line (default 0).
#' @param show_pre_treatment Logical: include pre-treatment periods (default TRUE).
#' @param facet_by_cohort Logical: facet by cohort in staggered mode (default FALSE).
#' @param color_significant Logical: color significant points differently (default TRUE).
#' @param ref_period Integer or NULL: reference period to normalize ATT to zero.
#' @param aggregation Character: "mean" or "weighted" aggregation method.
#' @param df_strategy Character: "conservative", "weighted", or "fallback".
#' @param return_data Logical: if TRUE, return plot data instead of ggplot object.
#' @param point_size Numeric: size of point markers (default 3).
#' @param errorbar_width Numeric: width of error bar caps (default 0.2).
#' @param title Character or NULL: custom plot title.
#' @param theme Function: ggplot2 theme function (default \code{ggplot2::theme_minimal}).
#' @param ... Additional arguments (currently unused).
#' @return A ggplot object, or a data.frame if \code{return_data = TRUE}.
#' @export
plot_event_study.lwdid_result <- function(x,
                                          ci_level = 0.95,
                                          reference_line = 0,
                                          show_pre_treatment = TRUE,
                                          facet_by_cohort = FALSE,
                                          color_significant = TRUE,
                                          ref_period = NULL,
                                          aggregation = c("mean", "weighted"),
                                          df_strategy = c("conservative",
                                                          "weighted",
                                                          "fallback"),
                                          return_data = FALSE,
                                          point_size = 3,
                                          errorbar_width = 0.2,
                                          title = NULL,
                                          theme = ggplot2::theme_minimal,
                                          ...) {

  # ==========================================================================
  # Step 1: Parameter validation
  # ==========================================================================

  # ggplot2 dependency check
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install via: install.packages('ggplot2')")
  }

  # Type checks
  if (!inherits(x, "lwdid_result")) {
    stop("x must be an lwdid_result object")
  }
  stopifnot(
    is.numeric(ci_level), length(ci_level) == 1L,
    ci_level > 0, ci_level < 1
  )
  if (!is.null(ref_period)) {
    stopifnot(
      is.numeric(ref_period), length(ref_period) == 1L,
      ref_period == as.integer(ref_period)
    )
    ref_period <- as.integer(ref_period)
  }
  aggregation <- match.arg(aggregation)
  df_strategy <- match.arg(df_strategy)
  if (!is.function(theme)) {
    stop("theme must be a ggplot2 theme function (e.g. theme_minimal)")
  }


  # ==========================================================================
  # Step 2: Data preparation — route by staggered vs common timing
  # ==========================================================================

  if (x$is_staggered) {
    # ------------------------------------------------------------------
    # Staggered mode
    # ------------------------------------------------------------------

    # Validate att_by_cohort_time
    cohort_df <- x$att_by_cohort_time
    if (is.null(cohort_df) || nrow(cohort_df) == 0L) {
      stop("att_by_cohort_time is empty, cannot create Event Study plot")
    }
    # Check required columns
    if (!all(c("att", "se") %in% names(cohort_df))) {
      missing_cols <- setdiff(c("att", "se"), names(cohort_df))
      stop(sprintf(
        "att_by_cohort_time missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ))
    }

    # event_time auto-calculation
    if (!"event_time" %in% names(cohort_df)) {
      if (!all(c("cohort", "period") %in% names(cohort_df))) {
        missing <- setdiff(c("cohort", "period"), names(cohort_df))
        stop(sprintf(
          paste0("att_by_cohort_time missing event_time column, ",
                 "and cannot compute automatically from period-cohort ",
                 "(missing columns: %s)"),
          paste(missing, collapse = ", ")
        ))
      }
      cohort_df$event_time <- cohort_df$period - cohort_df$cohort
    }
    cohort_df$`_source` <- "post_treatment"

    # Merge pre-treatment data BEFORE aggregation
    has_pre_data <- show_pre_treatment &&
      !is.null(x$att_pre_treatment) &&
      nrow(x$att_pre_treatment) > 0L
    if (!is.null(x$include_pretreatment) && !x$include_pretreatment) {
      has_pre_data <- FALSE
    }
    if (has_pre_data) {
      pre_df <- x$att_pre_treatment
      if (!"event_time" %in% names(pre_df) &&
          all(c("cohort", "period") %in% names(pre_df))) {
        pre_df$event_time <- pre_df$period - pre_df$cohort
      }
      pre_df$`_source` <- "pre_treatment"
      common_cols <- intersect(names(cohort_df), names(pre_df))
      cohort_df <- rbind(
        cohort_df[, common_cols, drop = FALSE],
        pre_df[, common_cols, drop = FALSE]
      )
    }

    # Filter out pre-treatment if not showing
    if (!show_pre_treatment) {
      cohort_df <- cohort_df[cohort_df$event_time >= 0L, , drop = FALSE]
    }


    # Validate cohort_sizes for weighted aggregation
    if (aggregation == "weighted") {
      if (is.null(x$cohort_sizes) || length(x$cohort_sizes) == 0L) {
        stop(paste0(
          "aggregation='weighted' requires x$cohort_sizes ",
          "(cohort sample sizes), but this field is empty"
        ))
      }
      # Check for non-positive values
      invalid <- x$cohort_sizes[x$cohort_sizes <= 0]
      if (length(invalid) > 0L) {
        invalid_str <- paste(
          sprintf("cohort %s=%s", names(invalid), invalid),
          collapse = ", "
        )
        stop(sprintf(
          paste0("cohort_sizes contains non-positive values: %s. ",
                 "All cohort sizes must be positive integers"),
          invalid_str
        ))
      }
    }

    # Aggregate by event time
    if (aggregation == "mean") {
      plot_data <- .aggregate_by_event_time(
        cohort_df, method = "mean", df_strategy = df_strategy
      )
      n_total_cohorts <- length(unique(cohort_df$cohort))
      if (n_total_cohorts > 1L) {
        warning(sprintf(
          paste0("Simple-average SE assumes %d cohorts are independent. ",
                 "When cohorts share control units, this assumption may ",
                 "lead to underestimated SE. Consider using ",
                 "aggregation='weighted' or checking inter-cohort independence."),
          n_total_cohorts
        ), call. = FALSE)
      }
    } else {
      plot_data <- .aggregate_by_event_time(
        cohort_df, method = "weighted",
        cohort_sizes = x$cohort_sizes,
        df_strategy = df_strategy
      )
    }

    # Sort by event_time
    plot_data <- plot_data[order(plot_data$event_time), , drop = FALSE]


  } else {
    # ------------------------------------------------------------------
    # Common Timing mode
    # ------------------------------------------------------------------
    plot_data <- x$att_by_period
    if (is.null(plot_data) || nrow(plot_data) == 0L) {
      stop("att_by_period is empty, cannot create Event Study plot")
    }
    # Check required columns
    if (!"att" %in% names(plot_data)) {
      stop("att_by_period missing required column: att")
    }
    if (!"se" %in% names(plot_data)) {
      stop("att_by_period missing required column: se")
    }

    # Add df_inference if missing
    if (!"df_inference" %in% names(plot_data)) {
      plot_data$df_inference <- x$df_inference
    }

    # event_time auto-calculation for Common Timing
    if (!"event_time" %in% names(plot_data) &&
        "period" %in% names(plot_data)) {
      t_time <- x$treatment_time %||% x$tpost1
      if (is.null(t_time)) {
        stop("Cannot compute event_time: missing treatment_time and tpost1")
      }
      plot_data$event_time <- plot_data$period - t_time
    }

    # Pre-treatment merge
    has_pre_data_ct <- show_pre_treatment &&
      !is.null(x$att_pre_treatment) &&
      nrow(x$att_pre_treatment) > 0L
    if (!is.null(x$include_pretreatment) && !x$include_pretreatment) {
      has_pre_data_ct <- FALSE
    }
    if (has_pre_data_ct) {
      pre_ct <- x$att_pre_treatment
      if (!"event_time" %in% names(pre_ct) &&
          "period" %in% names(pre_ct)) {
        t_time_pre <- x$treatment_time %||% x$tpost1
        pre_ct$event_time <- pre_ct$period - t_time_pre
      }
      if (!"df_inference" %in% names(pre_ct)) {
        pre_ct$df_inference <- x$df_inference
      }
      common_cols_ct <- intersect(names(plot_data), names(pre_ct))
      plot_data <- rbind(
        plot_data[, common_cols_ct, drop = FALSE],
        pre_ct[, common_cols_ct, drop = FALSE]
      )
    }

    # Filter pre-treatment
    if (!show_pre_treatment) {
      plot_data <- plot_data[plot_data$event_time >= 0L, , drop = FALSE]
    }

    # Add n_cohorts column for consistency
    if (!"n_cohorts" %in% names(plot_data)) {
      plot_data$n_cohorts <- NA_integer_
    }

    # Standardize columns to match staggered output format
    # (drop extra columns like period, t_stat, pvalue from att_by_period
    #  so that anchor row rbind in Step 4 works correctly)
    std_cols <- c("event_time", "att", "se", "df_inference", "n_cohorts")
    std_cols <- intersect(std_cols, names(plot_data))
    plot_data <- plot_data[, std_cols, drop = FALSE]

    # Sort
    plot_data <- plot_data[order(plot_data$event_time), , drop = FALSE]
  }


  # ==========================================================================
  # Step 3: Confidence interval calculation (per-event-time df)
  # ==========================================================================
  alpha_val <- 1 - ci_level
  t_crit <- stats::qt(1 - alpha_val / 2, df = plot_data$df_inference)
  plot_data$ci_lower <- plot_data$att - t_crit * plot_data$se
  plot_data$ci_upper <- plot_data$att + t_crit * plot_data$se

  # ==========================================================================
  # Step 4: Add anchor point
  # ==========================================================================
  anchor_time <- if (!is.null(ref_period)) ref_period else -1L

  # Pre-check: if show_pre_treatment=FALSE and ref_period is negative,
  # the ref_period would be in the filtered-out range. Error early.
  if (!is.null(ref_period) && !show_pre_treatment && ref_period < 0L) {
    available_et <- paste(
      sort(unique(plot_data$event_time)), collapse = ", "
    )
    stop(sprintf(
      paste0(
        "ref_period=%d is not in the data's event_time range ",
        "(pre-treatment data has been filtered out). ",
        "Set show_pre_treatment=TRUE ",
        "or choose a non-negative ref_period. ",
        "Available event_time: %s"
      ),
      ref_period, available_et
    ), call. = FALSE)
  }

  if (!(anchor_time %in% plot_data$event_time)) {
    # Anchor df: use minimum df from existing data, or Inf if empty
    anchor_df <- if (nrow(plot_data) > 0L) {
      min(plot_data$df_inference, na.rm = TRUE)
    } else {
      Inf
    }
    anchor_row <- data.frame(
      event_time   = anchor_time,
      att          = 0,
      se           = 0,
      df_inference = anchor_df,
      n_cohorts    = 0L,
      ci_lower     = 0,
      ci_upper     = 0,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, anchor_row)
    plot_data <- plot_data[order(plot_data$event_time), , drop = FALSE]
    rownames(plot_data) <- NULL
  }

  plot_data$is_anchor <- (plot_data$event_time == anchor_time)

  # ==========================================================================
  # Step 4b: Reference period normalization
  # ==========================================================================
  if (!is.null(ref_period)) {
    ref_row <- plot_data[plot_data$event_time == ref_period, , drop = FALSE]
    if (nrow(ref_row) == 0L) {
      # Check if filtered out by show_pre_treatment=FALSE
      is_filtered_out <- !show_pre_treatment && ref_period < 0L
      available_et <- paste(
        sort(unique(plot_data$event_time)), collapse = ", "
      )
      if (is_filtered_out) {
        stop(sprintf(
          paste0(
            "ref_period=%d is not in the data's event_time range ",
            "(pre-treatment data has been filtered out). ",
            "Set show_pre_treatment=TRUE ",
            "or choose a non-negative ref_period. ",
            "Available event_time: %s"
          ),
          ref_period, available_et
        ), call. = FALSE)
      } else {
        stop(sprintf(
          paste0(
            "ref_period=%d ",
            "is not in the data's event_time range. ",
            "Available event_time: %s"
          ),
          ref_period, available_et
        ), call. = FALSE)
      }
    }
    ref_att <- ref_row$att[1L]
    plot_data$att <- plot_data$att - ref_att
    plot_data$ci_lower <- plot_data$ci_lower - ref_att
    plot_data$ci_upper <- plot_data$ci_upper - ref_att
  }

  # ==========================================================================
  # Step 5: p-value, significance marking, and period type classification
  # ==========================================================================
  # Recompute p-value from t-distribution (R enhancement over Python)
  plot_data$pvalue <- 2 * stats::pt(
    -abs(plot_data$att / plot_data$se),
    df = plot_data$df_inference
  )
  # SE=0 or non-finite SE -> pvalue=NA (avoid division by zero)
  bad_se <- plot_data$se == 0 | !is.finite(plot_data$se)
  plot_data$pvalue[bad_se] <- NA_real_

  # Significance flag
  plot_data$significant <- !is.na(plot_data$pvalue) &
    plot_data$pvalue < (1 - ci_level)

  # Period type classification
  plot_data$period_type <- ifelse(
    plot_data$event_time < 0L,
    "pre_treatment",
    "post_treatment"
  )
  plot_data$period_type[plot_data$is_anchor] <- "anchor"


  # ==========================================================================
  # Step 6: Build ggplot layers
  # ==========================================================================

  # Base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$event_time, y = .data$att)) +
    ggplot2::geom_hline(yintercept = reference_line, linetype = "dashed", color = "gray50") +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      width = errorbar_width, color = "gray40"
    )

  if (color_significant) {
    # Mode A: significance color coding (R-specific)
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$significant, shape = .data$is_anchor),
        size = point_size
      ) +
      ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "#2166AC", "FALSE" = "gray60"),
        labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
        name = NULL
      )
  } else {
    # Mode B: pre/post treatment color (matches Python style)
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$period_type, shape = .data$is_anchor),
        size = point_size
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "pre_treatment"  = "gray60",
          "post_treatment" = "#2166AC",
          "anchor"         = "#E66100"
        ),
        name = NULL
      ) +
      ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none")
  }

  # Connection line
  p <- p + ggplot2::geom_line(alpha = 0.3, color = "gray40")

  # Pre/post treatment transition dotted line (excluding anchor points)
  pre_data <- plot_data[plot_data$event_time < 0L & !plot_data$is_anchor, , drop = FALSE]
  post_data <- plot_data[plot_data$event_time >= 0L & !plot_data$is_anchor, , drop = FALSE]
  if (nrow(pre_data) > 0L && nrow(post_data) > 0L) {
    last_pre <- pre_data[which.max(pre_data$event_time), , drop = FALSE]
    first_post <- post_data[which.min(post_data$event_time), , drop = FALSE]
    bridge_df <- data.frame(
      event_time = c(last_pre$event_time, first_post$event_time),
      att        = c(last_pre$att, first_post$att)
    )
    p <- p + ggplot2::geom_line(
      data = bridge_df,
      ggplot2::aes(x = .data$event_time, y = .data$att),
      linetype = "dotted", color = "gray50", linewidth = 0.8,
      inherit.aes = FALSE
    )
  }

  # Labels and subtitle
  effective_title <- if (!is.null(title)) title else "Event Study"
  subtitle_text <- .format_subtitle(
    aggregation, df_strategy, plot_data$df_inference, x$is_staggered
  )
  p <- p + ggplot2::labs(
    x        = "Event Time (Relative to Treatment)",
    y        = "ATT Estimate",
    title    = effective_title,
    subtitle = subtitle_text
  ) +
    theme()

  # x-axis integer ticks (one tick per event_time)
  p <- p + ggplot2::scale_x_continuous(
    breaks = sort(unique(plot_data$event_time))
  )

  # Remove minor grid lines
  p <- p + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  # Facet by cohort (staggered mode only, requires cohort column in plot_data)
  if (facet_by_cohort && "cohort" %in% names(plot_data)) {
    p <- p + ggplot2::facet_wrap(~ cohort)
  }


  # ==========================================================================
  # Step 7: Return result
  # ==========================================================================
  if (return_data) {
    return(list(plot = p, data = plot_data))
  } else {
    return(p)
  }
}
