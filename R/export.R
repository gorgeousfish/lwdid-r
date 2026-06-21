# ===========================================================================
# export.R — LaTeX and CSV export
#
# Implements result export framework: LaTeX table export (single model and
# multi-model comparison), CSV export, and LaTeX special character escaping.
# ===========================================================================

# ---------------------------------------------------------------------------
# .escape_latex: LaTeX special character escaping
# ---------------------------------------------------------------------------

#' Escape LaTeX special characters in text
#'
#' Internal helper that escapes LaTeX metacharacters such as backslashes,
#' underscores, hashes, dollar signs, percent signs, ampersands, braces,
#' tildes, and carets so that arbitrary text can be safely embedded in a
#' LaTeX document.
#'
#' @param text A single character string (or value coercible to character).
#'   \code{NULL} and \code{NA} are silently converted to \code{""}.
#' @return A single escaped character string.
#' @keywords internal
.escape_latex <- function(text) {
  # 1. NULL or NA → return ""
  if (is.null(text) || is.na(text)) return("")

  # 2. Non-character input → auto-convert
  if (!is.character(text)) text <- as.character(text)

  # 3. Backslash FIRST via placeholder (avoid double-escaping)
  text <- gsub("\\", "\x01BACKSLASH\x01", text, fixed = TRUE)

  # 4. Escape: _ # $ % & { }
  text <- gsub("([_#$%&{}])", "\\\\\\1", text)

  # 5. ~ → \textasciitilde{}
  text <- gsub("~", "\\\\textasciitilde{}", text)

  # 6. ^ → \textasciicircum{}
  text <- gsub("\\^", "\\\\textasciicircum{}", text)

  # 7. Replace placeholder with \textbackslash{}
  text <- gsub("\x01BACKSLASH\x01", "\\\\textbackslash{}", text)

  text
}

.export_estimator_description <- function(estimator) {
  if (is.null(estimator) || length(estimator) == 0L) {
    return("")
  }
  estimator <- as.character(estimator)[1L]
  if (is.na(estimator)) return("")
  key <- tolower(estimator)
  if (exists(".ESTIMATOR_DISPLAY_MAP", inherits = TRUE) &&
      key %in% names(.ESTIMATOR_DISPLAY_MAP)) {
    return(unname(.ESTIMATOR_DISPLAY_MAP[[key]]))
  }
  estimator
}

.export_vce_description <- function(vce_type, cluster_var = NULL) {
  if (is.null(vce_type) || length(vce_type) == 0L) {
    if (exists(".vce_description", inherits = TRUE)) {
      return(.vce_description(NULL, cluster_var))
    }
    return("")
  }
  vce_type <- as.character(vce_type)[1L]
  if (is.na(vce_type)) return("")
  if (exists(".vce_description", inherits = TRUE)) {
    return(.vce_description(vce_type, cluster_var))
  }
  vce_type
}

.wcb_export_metadata <- function(x) {
  details <- x$wcb_details
  if (is.null(details)) {
    return(NULL)
  }

  list(
    pvalue_wcb = x$pvalue_wcb %||% details$pvalue %||% NA_real_,
    se_wcb = x$se_wcb %||% details$se_bootstrap %||% NA_real_,
    ci_lower_wcb = x$ci_lower_wcb %||% details$ci_lower %||% NA_real_,
    ci_upper_wcb = x$ci_upper_wcb %||% details$ci_upper %||% NA_real_,
    wcb_weight_type = details$weight_type %||% NA_character_,
    wcb_n_clusters = details$n_clusters %||% NA_integer_,
    wcb_requested_n_bootstrap = details$requested_n_bootstrap %||% NA_integer_,
    wcb_actual_n_bootstrap = details$actual_n_bootstrap %||% NA_integer_,
    wcb_restricted_model = details$restricted_model %||% NA_character_,
    wcb_auto_triggered = isTRUE(x$wcb_auto_triggered)
  )
}

.format_latex_number <- function(x, digits) {
  if (length(x) == 0L || is.na(x[1L]) || !is.finite(x[1L])) {
    return("--")
  }
  sprintf("%.*f", digits, x[1L])
}

.format_latex_integer <- function(x) {
  if (length(x) == 0L || is.na(x[1L]) || !is.finite(x[1L])) {
    return("--")
  }
  sprintf("%d", as.integer(x[1L]))
}

.format_latex_pvalue <- function(x, digits) {
  if (length(x) == 0L || is.na(x[1L]) || !is.finite(x[1L])) {
    return("--")
  }
  if (x[1L] < 0.001) {
    return("$<$0.001")
  }
  sprintf("%.*f", digits, x[1L])
}

.format_latex_wcb_text <- function(x) {
  if (length(x) == 0L || is.null(x) || is.na(x[1L]) || !nzchar(as.character(x[1L]))) {
    return("--")
  }
  .escape_latex(as.character(x[1L]))
}

# ---------------------------------------------------------------------------
# S3 generics: to_latex / to_csv
# ---------------------------------------------------------------------------

#' Export an object to LaTeX table format
#'
#' S3 generic for converting estimation results to LaTeX table strings.
#'
#' @param x An object to export.
#' @param ... Additional arguments passed to methods.
#' @return A character string containing LaTeX table markup.
#' @seealso \code{\link{to_csv}}, \code{\link{to_dict}}
#' @export
to_latex <- function(x, ...) UseMethod("to_latex")

#' Export an object to CSV format
#'
#' S3 generic for writing estimation results to CSV files or returning
#' CSV-formatted strings.
#'
#' @param x An object to export.
#' @param ... Additional arguments passed to methods.
#' @return Depends on the method; typically a character string or
#'   invisible file path.
#' @seealso \code{\link{to_latex}}, \code{\link{to_dict}}
#' @export
to_csv <- function(x, ...) UseMethod("to_csv")

#' Convert an object to a dictionary (named list)
#'
#' @param x An object.
#' @param ... Additional arguments passed to methods.
#' @return A named list of parameters and statistics.
#' @seealso \code{\link{to_latex}}, \code{\link{to_csv}}
#' @export
to_dict <- function(x, ...) UseMethod("to_dict")

# ---------------------------------------------------------------------------
# to_latex.lwdid_result: LaTeX table export
# ---------------------------------------------------------------------------

#' @export
to_latex.lwdid_result <- function(x, caption = NULL, label = NULL,
                                    digits = 4L,
                                    include_ci = TRUE,
                                    include_diagnostics = FALSE,
                                    include_ri = TRUE,
                                    include_periods = FALSE,
                                    include_staggered = FALSE,
                                    include_event_time = FALSE,
                                    booktabs = TRUE,
                                    float = "htbp",
                                    se_in_parentheses = TRUE,
                                    stars = TRUE,
                                    file = NULL, ...) {
  # --- 0. Parameter validation ---
  stopifnot(is.numeric(digits), length(digits) == 1L, digits >= 0L)
  digits <- as.integer(digits)
  stopifnot(is.logical(include_ci), length(include_ci) == 1L)
  stopifnot(is.logical(include_diagnostics), length(include_diagnostics) == 1L)
  stopifnot(is.logical(include_ri), length(include_ri) == 1L)
  stopifnot(is.logical(include_periods), length(include_periods) == 1L)
  stopifnot(is.logical(include_staggered), length(include_staggered) == 1L)
  stopifnot(is.logical(include_event_time), length(include_event_time) == 1L)
  stopifnot(is.logical(booktabs), length(booktabs) == 1L)
  stopifnot(is.character(float), length(float) == 1L)
  stopifnot(is.logical(se_in_parentheses), length(se_in_parentheses) == 1L)
  stopifnot(is.logical(stars), length(stars) == 1L)
  if (!is.null(caption)) stopifnot(is.character(caption), length(caption) == 1L)
  if (!is.null(label)) stopifnot(is.character(label), length(label) == 1L)
  if (!is.null(file)) stopifnot(is.character(file), length(file) == 1L)

  # --- 1. Significance stars ---
  star_str <- ""
  if (stars && !is.na(x$pvalue)) {
    if (x$pvalue < 0.01) star_str <- "***"
    else if (x$pvalue < 0.05) star_str <- "**"
    else if (x$pvalue < 0.10) star_str <- "*"
  }

  rule_top <- if (booktabs) "\\toprule" else "\\hline"
  rule_mid <- if (booktabs) "\\midrule" else "\\hline"
  rule_bot <- if (booktabs) "\\bottomrule" else "\\hline"

  ci_level_pct <- round((1 - x$alpha) * 100)

  lines <- character(0)

  # --- 2. Table header ---
  lines <- c(lines, sprintf("\\begin{table}[%s]", float))
  lines <- c(lines, "\\centering")
  if (!is.null(caption)) {
    lines <- c(lines, sprintf("\\caption{%s}", .escape_latex(caption)))
  }
  if (!is.null(label)) {
    lines <- c(lines, sprintf("\\label{%s}", .escape_latex(label)))
  }
  lines <- c(lines, "\\begin{threeparttable}")
  lines <- c(lines, "\\begin{tabular}{lc}")
  lines <- c(lines, rule_top)

  # --- 3. Main body rows ---
  lines <- c(lines, sprintf("ATT & %s%s \\\\",
                              .format_latex_number(x$att, digits), star_str))

  if (se_in_parentheses) {
    lines <- c(lines, sprintf(" & (%s) \\\\",
                                .format_latex_number(x$se_att, digits)))
  } else {
    lines <- c(lines, sprintf("SE & %s \\\\",
                                .format_latex_number(x$se_att, digits)))
  }

  lines <- c(lines, rule_mid)

  if (include_ci) {
    lines <- c(lines, sprintf("CI [%d%%] & [%s, %s] \\\\",
                                ci_level_pct,
                                .format_latex_number(x$ci_lower, digits),
                                .format_latex_number(x$ci_upper, digits)))
  }

  # p-value row
  lines <- c(lines, sprintf("p-value & %s \\\\",
                              .format_latex_pvalue(x$pvalue, digits)))

  lines <- c(lines, rule_mid)

  # Sample info rows
  lines <- c(lines, sprintf("N & %s \\\\", .format_latex_integer(x$nobs)))
  lines <- c(lines, sprintf("$N_{treated}$ & %s \\\\",
                              .format_latex_integer(x$n_treated)))
  lines <- c(lines, sprintf("$N_{control}$ & %s \\\\",
                              .format_latex_integer(x$n_control)))
  lines <- c(lines, sprintf(
    "Estimator & %s \\\\",
    .escape_latex(.export_estimator_description(x$estimator))
  ))
  lines <- c(lines, sprintf(
    "VCE & %s \\\\",
    .escape_latex(.export_vce_description(x$vce_type, x$cluster_var))
  ))
  wcb_meta <- .wcb_export_metadata(x)
  if (!is.null(wcb_meta)) {
    lines <- c(lines, sprintf(
      "WCB weight & %s \\\\",
      .escape_latex(as.character(wcb_meta$wcb_weight_type))
    ))
    lines <- c(lines, sprintf(
      "WCB clusters & %s \\\\",
      .format_latex_integer(wcb_meta$wcb_n_clusters)
    ))
    lines <- c(lines, sprintf(
      "WCB draws & %s actual / %s requested \\\\",
      .format_latex_integer(wcb_meta$wcb_actual_n_bootstrap),
      .format_latex_integer(wcb_meta$wcb_requested_n_bootstrap)
    ))
    if (!is.na(wcb_meta$wcb_restricted_model) &&
        nzchar(wcb_meta$wcb_restricted_model)) {
      lines <- c(lines, sprintf(
        "WCB restricted model & %s \\\\",
        .escape_latex(as.character(wcb_meta$wcb_restricted_model))
      ))
    }
  }
  lines <- c(lines, sprintf("Rolling & %s \\\\",
                              .escape_latex(as.character(x$rolling))))

  # --- 4. Diagnostics (optional) ---
  if (include_diagnostics && !is.null(x$diagnostics)) {
    lines <- c(lines, rule_mid)
    if (!is.null(x$diagnostics$trends)) {
      lines <- c(lines, sprintf("Recommendation & %s \\\\",
                                  .escape_latex(x$diagnostics$trends$recommended_method)))
    }
    if (!is.null(x$diagnostics$selection)) {
      lines <- c(lines, sprintf("Selection Risk & %s \\\\",
                                  .escape_latex(x$diagnostics$selection$selection_risk)))
    }
  }

  # --- 5. RI info (optional) ---
  if (include_ri && !is.null(x$ri_pvalue)) {
    lines <- c(lines, rule_mid)
    lines <- c(lines, sprintf("RI p-value & %s \\\\",
                                .format_latex_pvalue(x$ri_pvalue, digits)))
    if (!is.null(x$ri_observed_stat)) {
      lines <- c(lines, sprintf("RI observed statistic & %s \\\\",
                                  .format_latex_number(x$ri_observed_stat, digits)))
    }
    if (!is.null(x$ri_estimator)) {
      lines <- c(lines, sprintf("RI estimator & %s \\\\",
                                  .format_latex_wcb_text(x$ri_estimator)))
    }
    if (!is.null(x$ri_seed)) {
      lines <- c(lines, sprintf("RI seed & %s \\\\",
                                  .format_latex_integer(x$ri_seed)))
    }
    if (!is.null(x$rireps)) {
      lines <- c(lines, sprintf("RI reps & %s \\\\",
                                  .format_latex_integer(x$rireps)))
    }
  }

  # --- 6. Main table footer ---
  lines <- c(lines, rule_bot)
  lines <- c(lines, "\\end{tabular}")
  if (stars) {
    lines <- c(lines, "\\begin{tablenotes}\\small")
    lines <- c(lines, "\\item \\textit{Note:} $^{***}$p$<$0.01, $^{**}$p$<$0.05, $^{*}$p$<$0.10")
    lines <- c(lines, "\\end{tablenotes}")
  }
  lines <- c(lines, "\\end{threeparttable}")

  # --- 7. Period effects table (optional) ---
  if (include_periods && !is.null(x$att_by_period) &&
      is.data.frame(x$att_by_period) && nrow(x$att_by_period) > 0) {
    period_rows <- .lwdid_effect_rows_with_supplied_vcov(
      x$att_by_period, x, "by_period"
    )
    lines <- c(lines, "\\vspace{1em}")
    lines <- c(lines, "\\begin{threeparttable}")
    lines <- c(lines, "\\begin{tabular}{lccccc}")
    lines <- c(lines, rule_top)
    lines <- c(lines, "Period & ATT & SE & t-stat & p-value & CI \\\\")
    lines <- c(lines, rule_mid)

    for (i in seq_len(nrow(period_rows))) {
      row <- period_rows[i, ]

      # Period identifier (period or event_time)
      period_id <- if ("period" %in% names(row)) row$period
                   else if ("event_time" %in% names(row)) row$event_time
                   else i
      period_label <- .escape_latex(as.character(period_id))

      # ATT and SE
      att_val <- .format_latex_number(row$att, digits)
      se_val <- .format_latex_number(row$se, digits)

      # Period stars
      period_star <- ""
      if (stars && "pvalue" %in% names(row) && !is.na(row$pvalue)) {
        if (row$pvalue < 0.01) period_star <- "***"
        else if (row$pvalue < 0.05) period_star <- "**"
        else if (row$pvalue < 0.10) period_star <- "*"
      }

      # t-statistic
      t_val <- if ("t_stat" %in% names(row)) {
        .format_latex_number(row$t_stat, digits)
      } else if ("t" %in% names(row)) {
        .format_latex_number(row$t, digits)
      } else {
        .format_latex_number(row$att / row$se, digits)
      }

      # p-value
      p_str <- if (!"pvalue" %in% names(row)) "--"
               else .format_latex_pvalue(row$pvalue, digits)

      # CI (handle column name variants)
      ci_lo <- if ("ci_lower" %in% names(row)) row$ci_lower
               else if ("conf.low" %in% names(row)) row[["conf.low"]]
               else NA
      ci_hi <- if ("ci_upper" %in% names(row)) row$ci_upper
               else if ("conf.high" %in% names(row)) row[["conf.high"]]
               else NA
      ci_str <- if (is.na(ci_lo) || is.na(ci_hi) ||
                    !is.finite(ci_lo) || !is.finite(ci_hi)) "--"
                else sprintf("[%s, %s]",
                             .format_latex_number(ci_lo, digits),
                             .format_latex_number(ci_hi, digits))

      lines <- c(lines, sprintf("%s & %s%s & %s & %s & %s & %s \\\\",
                                  period_label, att_val, period_star,
                                  se_val, t_val, p_str, ci_str))
    }

    lines <- c(lines, rule_bot)
    lines <- c(lines, "\\end{tabular}")
    lines <- c(lines, "\\end{threeparttable}")
  }

  # --- 8. Staggered cohort effects table (optional) ---
  if (include_staggered && isTRUE(x$is_staggered) &&
      !is.null(x$att_by_cohort) && is.data.frame(x$att_by_cohort) &&
      nrow(x$att_by_cohort) > 0) {
    cohort_rows <- .lwdid_effect_rows_with_supplied_vcov(
      x$att_by_cohort, x, "by_cohort"
    )
    has_effective_weight_note <- FALSE
    if (!is.null(x$cohort_weights) && !is.null(x$effective_weights)) {
      shared_weight_ids <- intersect(names(x$cohort_weights), names(x$effective_weights))
      if (length(shared_weight_ids) > 0L) {
        computation_weight_vals <- as.numeric(x$cohort_weights[shared_weight_ids])
        effective_weight_vals <- as.numeric(x$effective_weights[shared_weight_ids])
        valid_weight_vals <- is.finite(computation_weight_vals) &
          is.finite(effective_weight_vals)
        has_effective_weight_note <- any(
          abs(computation_weight_vals[valid_weight_vals] -
                effective_weight_vals[valid_weight_vals]) > 1e-12
        )
      }
    }

    lines <- c(lines, "\\vspace{1em}")
    lines <- c(lines, "\\begin{threeparttable}")
    lines <- c(lines, "\\begin{tabular}{lcccc}")
    lines <- c(lines, rule_top)
    lines <- c(lines, "Cohort & ATT & SE & N & Weight \\\\")
    lines <- c(lines, rule_mid)

    for (i in seq_len(nrow(cohort_rows))) {
      row <- cohort_rows[i, ]
      cohort_id <- if ("cohort" %in% names(row)) row$cohort else i
      cohort_label <- .escape_latex(as.character(cohort_id))
      cohort_att <- .format_latex_number(row$att, digits)
      cohort_se <- .format_latex_number(row$se, digits)
      cohort_n <- if ("n" %in% names(row)) .format_latex_integer(row$n)
                  else if ("nobs" %in% names(row)) .format_latex_integer(row$nobs)
                  else {
                    cohort_sizes <- x$cohort_sample_sizes %||% x$cohort_sizes
                    cid <- as.character(cohort_id)
                    if (!is.null(cohort_sizes) && cid %in% names(cohort_sizes)) {
                      .format_latex_integer(cohort_sizes[[cid]])
                    } else {
                      "--"
                    }
                  }

      # Weight: prefer att_by_cohort$weight, fallback to cohort_weights
      cohort_w <- if ("weight" %in% names(row)) {
        .format_latex_number(row$weight, digits)
      } else if (!is.null(x$cohort_weights)) {
        cid <- as.character(cohort_id)
        if (cid %in% names(x$cohort_weights)) {
          .format_latex_number(x$cohort_weights[[cid]], digits)
        } else "--"
      } else "--"

      # Cohort stars (if pvalue column exists)
      cohort_star <- ""
      if (stars && "pvalue" %in% names(row) && !is.na(row$pvalue)) {
        if (row$pvalue < 0.01) cohort_star <- "***"
        else if (row$pvalue < 0.05) cohort_star <- "**"
        else if (row$pvalue < 0.10) cohort_star <- "*"
      }

      lines <- c(lines, sprintf("%s & %s%s & %s & %s & %s \\\\",
                                  cohort_label, cohort_att, cohort_star,
                                  cohort_se, cohort_n, cohort_w))
    }

    lines <- c(lines, rule_mid)
    lines <- c(lines, sprintf("Overall (weighted) & %s%s & %s & %s & -- \\\\",
                                .format_latex_number(x$att, digits), star_str,
                                .format_latex_number(x$se_att, digits),
                                .format_latex_integer(x$nobs)))
    lines <- c(lines, rule_bot)
    lines <- c(lines, "\\end{tabular}")
    if (has_effective_weight_note) {
      lines <- c(lines, "\\begin{tablenotes}\\footnotesize")
      lines <- c(lines, paste0(
        "\\item \\textit{Note:} Weight reports the aggregation weight used ",
        "for the overall effect. Effective post-dropna regression-sample ",
        "weights are available in CSV and dictionary exports."
      ))
      lines <- c(lines, "\\end{tablenotes}")
    }
    lines <- c(lines, "\\end{threeparttable}")
  }

  # --- 9. Event-time effects table (optional) ---
  if (include_event_time && isTRUE(x$is_staggered) &&
      identical(x$aggregate, "event_time") &&
      !is.null(x$event_time_effects) && length(x$event_time_effects) > 0L) {
    event_rows <- extract_effects(x, type = "event_time")
    if (nrow(event_rows) > 0L) {
      has_support_metadata <- all(c("min_n_treated", "max_weight_cv") %in%
                                    names(event_rows))
      se_aggregation <- if ("se_aggregation" %in% names(event_rows)) {
        unique(stats::na.omit(event_rows$se_aggregation))
      } else {
        character(0)
      }
      covariance_assumption <- if ("covariance_assumption" %in% names(event_rows)) {
        unique(stats::na.omit(event_rows$covariance_assumption))
      } else {
        character(0)
      }

      lines <- c(lines, "\\vspace{1em}")
      lines <- c(lines, "\\begin{threeparttable}")
      lines <- c(lines, "\\begin{tabular}{lccccccc}")
      lines <- c(lines, rule_top)
      lines <- c(
        lines,
        "Event time & WATT(e) & SE & p-value & CI & Cohorts & Min treated & Max weight CV \\\\"
      )
      lines <- c(lines, rule_mid)

      for (i in seq_len(nrow(event_rows))) {
        row <- event_rows[i, , drop = FALSE]
        event_star <- ""
        if (stars && "pvalue" %in% names(row) && !is.na(row$pvalue)) {
          if (row$pvalue < 0.01) event_star <- "***"
          else if (row$pvalue < 0.05) event_star <- "**"
          else if (row$pvalue < 0.10) event_star <- "*"
        }
        ci_str <- if (all(c("ci_lower", "ci_upper") %in% names(row)) &&
                      is.finite(row$ci_lower) && is.finite(row$ci_upper)) {
          sprintf(
            "[%s, %s]",
            .format_latex_number(row$ci_lower, digits),
            .format_latex_number(row$ci_upper, digits)
          )
        } else {
          "--"
        }
        n_cohorts <- if ("n_cohorts" %in% names(row)) {
          .format_latex_integer(row$n_cohorts)
        } else {
          "--"
        }
        min_treated <- if ("min_n_treated" %in% names(row)) {
          .format_latex_integer(row$min_n_treated)
        } else {
          "--"
        }
        max_weight_cv <- if ("max_weight_cv" %in% names(row)) {
          .format_latex_number(row$max_weight_cv, digits)
        } else {
          "--"
        }
        p_str <- if ("pvalue" %in% names(row)) {
          .format_latex_pvalue(row$pvalue, digits)
        } else {
          "--"
        }
        lines <- c(lines, sprintf(
          "%s & %s%s & %s & %s & %s & %s & %s & %s \\\\",
          .escape_latex(as.character(row$event_time)),
          .format_latex_number(row$att, digits),
          event_star,
          .format_latex_number(row$se, digits),
          p_str,
          ci_str,
          n_cohorts,
          min_treated,
          max_weight_cv
        ))
      }

      lines <- c(lines, rule_bot)
      lines <- c(lines, "\\end{tabular}")
      lines <- c(lines, "\\begin{tablenotes}\\footnotesize")
      note <- paste0(
        "\\item \\textit{Note:} WATT(e) rows are event-time weighted ",
        "averages of available cohort-period effects."
      )
      if (length(se_aggregation) == 1L &&
          identical(se_aggregation, "diagonal_weighted_cohort_se")) {
        note <- paste0(
          note,
          " Standard errors use diagonal cohort-SE aggregation."
        )
      }
      if (length(covariance_assumption) == 1L &&
          identical(covariance_assumption, "zero_cross_cohort_covariance")) {
        note <- paste0(note, " Cross-cohort covariance is not modeled.")
      }
      if (has_support_metadata) {
        note <- paste0(
          note,
          " Min treated and max weight CV summarize contributing cells."
        )
      }
      lines <- c(lines, note)
      lines <- c(lines, "\\end{tablenotes}")
      lines <- c(lines, "\\end{threeparttable}")
    }
  }

  # --- 9. Table end and output ---
  lines <- c(lines, "\\end{table}")
  result <- paste(lines, collapse = "\n")

  if (!is.null(file)) {
    writeLines(result, file)
    message(sprintf("LaTeX table written to: %s", file))
    return(invisible(result))
  }

  result
}

# ---------------------------------------------------------------------------
# to_csv.lwdid_result: CSV export
# ---------------------------------------------------------------------------

.enrich_cohort_csv_metadata <- function(df, x) {
  if (!is.data.frame(df) || !"cohort" %in% names(df)) {
    return(df)
  }

  cohort_ids <- as.character(df$cohort)
  if (!"n" %in% names(df) && !"nobs" %in% names(df)) {
    cohort_sizes <- x$cohort_sample_sizes %||% x$cohort_sizes
    if (!is.null(cohort_sizes)) {
      df$n <- as.integer(cohort_sizes[cohort_ids])
    }
  }
  if (!"weight" %in% names(df) && !is.null(x$cohort_weights)) {
    df$weight <- as.numeric(x$cohort_weights[cohort_ids])
  }
  if (!"effective_weight" %in% names(df) && !is.null(x$effective_weights)) {
    df$effective_weight <- as.numeric(x$effective_weights[cohort_ids])
  }
  df
}

.format_export_reason_counts <- function(reason_counts) {
  if (is.null(reason_counts) ||
      !is.data.frame(reason_counts) ||
      nrow(reason_counts) == 0L) {
    return("")
  }
  paste(
    sprintf("%s: %d", reason_counts$reason, as.integer(reason_counts$n)),
    collapse = ", "
  )
}

#' @export
to_csv.lwdid_result <- function(x, file, what = "summary", ...) {
  # --- 0. Parameter validation ---
  stopifnot(is.character(file), length(file) == 1L)
  what <- match.arg(
    what,
    c(
      "summary", "by_period", "by_cohort", "all",
      "event_time_contributions"
    )
  )

  # --- 1. Build data.frame based on what ---
  df <- switch(what,
    "summary" = {
      skipped_diagnostics <- .validate_skipped_diagnostics(
        x$skipped_summary,
        x$skipped_pairs
      )
      event_time_omissions <- .validate_event_time_omissions(
        x$event_time_omission_summary,
        x$event_time_omissions
      )
      base_df <- data.frame(
        att = x$att,
        se = x$se_att,
        t_stat = x$t_stat,
        pvalue = x$pvalue,
        ci_lower = x$ci_lower,
        ci_upper = x$ci_upper,
        alpha = x$alpha,
        nobs = x$nobs,
        n_treated = x$n_treated,
        n_control = x$n_control,
        estimator = x$estimator %||% NA_character_,
        vce_type = x$vce_type %||% NA_character_,
        cluster_var = x$cluster_var %||% NA_character_,
        n_clusters = x$n_clusters %||% NA_integer_,
        rolling = as.character(x$rolling %||% NA_character_),
        method = x$method %||% NA_character_,
        is_staggered = isTRUE(x$is_staggered),
        aggregate = x$aggregate %||% NA_character_,
        control_group_used = x$control_group_used %||% NA_character_,
        n_cohorts = x$n_cohorts %||% NA_integer_,
        n_never_treated = x$n_never_treated %||% NA_integer_,
        n_skipped_pairs = skipped_diagnostics$n_skipped_pairs,
        skipped_summary = .format_export_reason_counts(
          skipped_diagnostics$skipped_summary
        ),
        n_omitted_event_times = event_time_omissions$n_omitted_event_times,
        event_time_omission_summary = .format_export_reason_counts(
          event_time_omissions$event_time_omission_summary
        ),
        stringsAsFactors = FALSE
      )
      if (!is.null(x$ri_pvalue)) {
        base_df$ri_pvalue <- x$ri_pvalue
        base_df$ri_observed_stat <- x$ri_observed_stat %||% NA_real_
        base_df$ri_estimator <- x$ri_estimator %||% NA_character_
        base_df$ri_seed <- if (!is.null(x$ri_seed)) x$ri_seed else NA
        base_df$rireps <- if (!is.null(x$rireps)) x$rireps else NA
      }
      event_time_support <- .summarize_event_time_support_metadata(x)
      if (!is.null(event_time_support)) {
        for (nm in names(event_time_support)) {
          base_df[[paste0("event_time_", nm)]] <- event_time_support[[nm]]
        }
      }
      event_time_bootstrap <- .summarize_event_time_bootstrap_metadata(x)
      if (!is.null(event_time_bootstrap)) {
        for (nm in names(event_time_bootstrap)) {
          base_df[[paste0("event_time_bootstrap_", nm)]] <-
            event_time_bootstrap[[nm]]
        }
      }
      wcb_meta <- .wcb_export_metadata(x)
      if (!is.null(wcb_meta)) {
        for (nm in names(wcb_meta)) {
          base_df[[nm]] <- wcb_meta[[nm]]
        }
      }
      base_df
    },
    "by_period" = {
      if (is.null(x$att_by_period) ||
          (is.data.frame(x$att_by_period) && nrow(x$att_by_period) == 0L)) {
        stop("No period-specific results available")
      }
      .lwdid_effect_rows_with_supplied_vcov(x$att_by_period, x, "by_period")
    },
    "by_cohort" = {
      if (is.null(x$att_by_cohort) ||
          (is.data.frame(x$att_by_cohort) && nrow(x$att_by_cohort) == 0L)) {
        stop("No cohort-specific results (Staggered mode only).",
             call. = FALSE)
      }
      cohort_rows <- .lwdid_effect_rows_with_supplied_vcov(
        x$att_by_cohort, x, "by_cohort"
      )
      .enrich_cohort_csv_metadata(cohort_rows, x)
    },
    "all" = {
      if (!is.null(x$att_by_cohort_time)) {
        .enrich_cohort_csv_metadata(x$att_by_cohort_time, x)
      }
      else if (!is.null(x$att_by_period)) {
        .lwdid_effect_rows_with_supplied_vcov(
          x$att_by_period, x, "by_period"
        )
      }
      else data.frame(att = x$att, se = x$se_att)
    },
    "event_time_contributions" = {
      extract_effects(x, type = "event_time_contributions")
    }
  )

  # --- 2. Write CSV ---
  utils::write.csv(df, file = file, row.names = FALSE, ...)
  message(sprintf("CSV written to: %s (%d rows)", file, nrow(df)))
  invisible(df)
}

# ---------------------------------------------------------------------------
# to_dict.lwdid_result: dictionary export
# ---------------------------------------------------------------------------

#' @export
to_dict.lwdid_result <- function(x, ...) {
  warnings_count <- if (is.null(x$warnings_log)) {
    0L
  } else if (is.data.frame(x$warnings_log)) {
    nrow(x$warnings_log)
  } else {
    length(x$warnings_log)
  }
  cohort_sample_sizes <- x$cohort_sample_sizes %||% x$cohort_sizes
  call_val <- if (!is.null(x$call)) {
    x$call
  } else if (!is.null(x$metadata) && !is.null(x$metadata$call)) {
    x$metadata$call
  } else {
    NULL
  }
  skipped_diagnostics <- .validate_skipped_diagnostics(
    x$skipped_summary,
    x$skipped_pairs
  )
  event_time_omissions <- .validate_event_time_omissions(
    x$event_time_omission_summary,
    x$event_time_omissions
  )
  event_time_support <- .summarize_event_time_support_metadata(x)
  event_time_bootstrap <- .summarize_event_time_bootstrap_metadata(x)

  out <- list(
    att = x$att,
    se_att = x$se_att,
    t_stat = x$t_stat,
    pvalue = x$pvalue,
    ci_lower = x$ci_lower,
    ci_upper = x$ci_upper,
    df_resid = x$df_resid,
    df_inference = x$df_inference,
    nobs = x$nobs,
    n_treated = x$n_treated,
    n_control = x$n_control,
    K = x$K,
    tpost1 = x$tpost1,
    depvar = x$depvar,
    rolling = x$rolling,
    vce_type = x$vce_type,
    cluster_var = x$cluster_var,
    n_clusters = x$n_clusters,
    estimator = x$estimator,
    method = x$method,
    alpha = x$alpha,
    is_staggered = x$is_staggered,
    controls_used = x$controls_used,
    controls = x$controls,
    include_pretreatment = x$include_pretreatment,
    control_group = x$control_group,
    control_group_used = x$control_group_used,
    aggregate = x$aggregate,
    cohorts = x$cohorts,
    cohort_sizes = x$cohort_sizes,
    cohort_sample_sizes = cohort_sample_sizes,
    cohort_weights = x$cohort_weights,
    effective_weights = x$effective_weights,
    n_never_treated = x$n_never_treated,
    n_units = x$n_units,
    n_periods = x$n_periods,
    n_cohorts = x$n_cohorts,
    skipped_pairs = skipped_diagnostics$skipped_pairs,
    skipped_summary = skipped_diagnostics$skipped_summary,
    n_skipped_pairs = skipped_diagnostics$n_skipped_pairs,
    event_time_omissions = event_time_omissions$event_time_omissions,
    event_time_omission_summary =
      event_time_omissions$event_time_omission_summary,
    n_omitted_event_times = event_time_omissions$n_omitted_event_times,
    event_time_support_summary = event_time_support,
    event_time_bootstrap_summary = event_time_bootstrap,
    att_overall = x$att_overall,
    se_overall = x$se_overall,
    ri_pvalue = x$ri_pvalue,
    ri_seed = x$ri_seed,
    rireps = x$rireps,
    ri_method = x$ri_method,
    ri_valid = x$ri_valid,
    ri_failed = x$ri_failed,
    ri_error = x$ri_error,
    ri_target = x$ri_target,
    ri_observed_stat = x$ri_observed_stat,
    ri_estimator = x$ri_estimator,
    warnings_count = warnings_count,
    call = call_val,
    lwdid_version = x$lwdid_version,
    ivar = x$ivar,
    tvar = x$tvar,
    is_quarterly = x$is_quarterly
  )
  wcb_meta <- .wcb_export_metadata(x)
  if (!is.null(wcb_meta)) {
    out <- c(out, wcb_meta)
  }
  if (!is.null(x$event_time_effects) &&
      !is.data.frame(x$event_time_effects) &&
      length(x$event_time_effects) > 0L) {
    event_time_contributions <- extract_effects(
      x,
      type = "event_time_contributions"
    )
    if (nrow(event_time_contributions) > 0L) {
      out$event_time_contributions <- event_time_contributions
    }
  }
  out
}

# ---------------------------------------------------------------------------
# to_latex_comparison: Multi-model comparison LaTeX table
# ---------------------------------------------------------------------------

#' Export multiple lwdid results to a comparison LaTeX table
#'
#' @param ... One or more objects of class \code{lwdid_result}.
#' @param model_names Optional character vector of column labels. Must have the
#'   same length as the number of supplied models.
#' @param caption Optional table caption.
#' @param label Optional LaTeX label.
#' @param digits Non-negative integer number of digits to print.
#' @param include_ci Logical; if \code{TRUE}, include confidence intervals.
#' @param booktabs Logical; if \code{TRUE}, use \code{booktabs} rules.
#' @param stars Logical; if \code{TRUE}, append significance stars.
#' @param file Optional output path. When supplied, the rendered LaTeX string
#'   is also written to disk.
#' @return A single character string containing LaTeX table markup. When
#'   \code{file} is supplied, the same string is returned invisibly after
#'   writing the file.
#' @export
to_latex_comparison <- function(..., model_names = NULL,
                                 caption = NULL, label = NULL,
                                 digits = 4L, include_ci = TRUE,
                                 booktabs = TRUE, stars = TRUE,
                                 file = NULL) {
  # --- 0. Parameter validation ---
  models <- list(...)
  n_models <- length(models)
  if (n_models == 0L) stop("At least one model is required")

  stopifnot(is.numeric(digits), length(digits) == 1L, digits >= 0L)
  digits <- as.integer(digits)
  stopifnot(is.logical(include_ci), length(include_ci) == 1L)
  stopifnot(is.logical(booktabs), length(booktabs) == 1L)
  stopifnot(is.logical(stars), length(stars) == 1L)
  if (!is.null(model_names)) {
    stopifnot(is.character(model_names), length(model_names) == n_models)
  }
  if (!is.null(caption)) stopifnot(is.character(caption), length(caption) == 1L)
  if (!is.null(label)) stopifnot(is.character(label), length(label) == 1L)
  if (!is.null(file)) stopifnot(is.character(file), length(file) == 1L)

  # Validate all models are lwdid_result
  for (m in models) {
    if (!inherits(m, "lwdid_result")) {
      stop("All models must be lwdid_result objects")
    }
  }

  # --- 1. Model names ---
  if (is.null(model_names)) {
    nm <- names(models)
    if (!is.null(nm) && all(nzchar(nm))) {
      model_names <- nm
    } else {
      model_names <- paste0("(", seq_len(n_models), ")")
    }
  }

  # --- 2. Significance stars ---
  if (stars) {
    star_strs <- vapply(models, function(m) {
      if (is.na(m$pvalue)) return("")
      if (m$pvalue < 0.01) return("***")
      if (m$pvalue < 0.05) return("**")
      if (m$pvalue < 0.10) return("*")
      ""
    }, character(1))
  } else {
    star_strs <- rep("", n_models)
  }

  if (include_ci) {
    alpha_vals <- vapply(models, function(m) m$alpha, numeric(1))
    if (any(!is.finite(alpha_vals))) {
      stop("All models must have finite alpha values when include_ci = TRUE")
    }
    if (length(unique(alpha_vals)) != 1L) {
      stop("All models must use the same alpha when include_ci = TRUE")
    }
  }

  rule_top <- if (booktabs) "\\toprule" else "\\hline"
  rule_mid <- if (booktabs) "\\midrule" else "\\hline"
  rule_bot <- if (booktabs) "\\bottomrule" else "\\hline"

  col_spec <- paste0("l", paste(rep("c", n_models), collapse = ""))

  lines <- character(0)

  # --- 3. Table header ---
  lines <- c(lines, "\\begin{table}[htbp]")
  lines <- c(lines, "\\centering")
  if (!is.null(caption)) {
    lines <- c(lines, sprintf("\\caption{%s}", .escape_latex(caption)))
  }
  if (!is.null(label)) {
    lines <- c(lines, sprintf("\\label{%s}", .escape_latex(label)))
  }
  if (stars) lines <- c(lines, "\\begin{threeparttable}")
  lines <- c(lines, sprintf("\\begin{tabular}{%s}", col_spec))
  lines <- c(lines, rule_top)

  # --- 4. Header row ---
  escaped_names <- vapply(model_names, .escape_latex, character(1))
  header_cells <- paste(escaped_names, collapse = " & ")
  lines <- c(lines, sprintf(" & %s \\\\", header_cells))
  lines <- c(lines, rule_mid)

  # --- 5. ATT row ---
  att_cells <- vapply(seq_len(n_models), function(i) {
    sprintf("%s%s", .format_latex_number(models[[i]]$att, digits),
            star_strs[i])
  }, character(1))
  lines <- c(lines, sprintf("ATT & %s \\\\", paste(att_cells, collapse = " & ")))

  # --- 6. SE row (always parentheses) ---
  se_cells <- vapply(models, function(m) {
    sprintf("(%s)", .format_latex_number(m$se_att, digits))
  }, character(1))
  lines <- c(lines, sprintf(" & %s \\\\", paste(se_cells, collapse = " & ")))

  # --- 7. CI row (optional) ---
  if (include_ci) {
    ci_cells <- vapply(models, function(m) {
      ci_level_pct <- round((1 - m$alpha) * 100)
      sprintf("[%s, %s]",
              .format_latex_number(m$ci_lower, digits),
              .format_latex_number(m$ci_upper, digits))
    }, character(1))
    ci_level_pct <- round((1 - models[[1]]$alpha) * 100)
    lines <- c(lines, sprintf("CI [%d%%] & %s \\\\",
                                ci_level_pct,
                                paste(ci_cells, collapse = " & ")))
  }

  # --- 8. p-value row ---
  pval_cells <- vapply(models, function(m) {
    .format_latex_pvalue(m$pvalue, digits)
  }, character(1))
  lines <- c(lines, sprintf("p-value & %s \\\\",
                              paste(pval_cells, collapse = " & ")))
  lines <- c(lines, rule_mid)

  # --- 9. Sample info rows ---
  n_cells <- vapply(models, function(m) .format_latex_integer(m$nobs),
                    character(1))
  lines <- c(lines, sprintf("N & %s \\\\", paste(n_cells, collapse = " & ")))

  est_cells <- vapply(models, function(m) {
    .escape_latex(.export_estimator_description(m$estimator))
  }, character(1))
  lines <- c(lines, sprintf("Estimator & %s \\\\",
                              paste(est_cells, collapse = " & ")))

  vce_cells <- vapply(models, function(m) {
    .escape_latex(.export_vce_description(m$vce_type, m$cluster_var))
  }, character(1))
  lines <- c(lines, sprintf("VCE & %s \\\\",
                              paste(vce_cells, collapse = " & ")))

  wcb_meta <- lapply(models, .wcb_export_metadata)
  has_wcb <- any(vapply(wcb_meta, Negate(is.null), logical(1)))
  if (has_wcb) {
    wcb_cell <- function(meta, field, formatter) {
      if (is.null(meta)) {
        return("--")
      }
      formatter(meta[[field]])
    }

    wcb_weight_cells <- vapply(
      wcb_meta,
      wcb_cell,
      character(1),
      field = "wcb_weight_type",
      formatter = .format_latex_wcb_text
    )
    lines <- c(lines, sprintf("WCB weight & %s \\\\",
                              paste(wcb_weight_cells, collapse = " & ")))

    wcb_cluster_cells <- vapply(
      wcb_meta,
      wcb_cell,
      character(1),
      field = "wcb_n_clusters",
      formatter = .format_latex_integer
    )
    lines <- c(lines, sprintf("WCB clusters & %s \\\\",
                              paste(wcb_cluster_cells, collapse = " & ")))

    wcb_draw_cells <- vapply(wcb_meta, function(meta) {
      if (is.null(meta)) {
        return("--")
      }
      sprintf(
        "%s actual / %s requested",
        .format_latex_integer(meta$wcb_actual_n_bootstrap),
        .format_latex_integer(meta$wcb_requested_n_bootstrap)
      )
    }, character(1))
    lines <- c(lines, sprintf("WCB draws & %s \\\\",
                              paste(wcb_draw_cells, collapse = " & ")))
  }

  roll_cells <- vapply(models, function(m) {
    .escape_latex(as.character(m$rolling))
  }, character(1))
  lines <- c(lines, sprintf("Rolling & %s \\\\",
                              paste(roll_cells, collapse = " & ")))

  # --- 10. Table footer ---
  lines <- c(lines, rule_bot)
  lines <- c(lines, "\\end{tabular}")
  if (stars) {
    lines <- c(lines, "\\begin{tablenotes}\\small")
    lines <- c(lines, "\\item \\textit{Note:} $^{***}$p$<$0.01, $^{**}$p$<$0.05, $^{*}$p$<$0.10")
    lines <- c(lines, "\\end{tablenotes}")
    lines <- c(lines, "\\end{threeparttable}")
  }
  lines <- c(lines, "\\end{table}")

  result <- paste(lines, collapse = "\n")

  if (!is.null(file)) {
    writeLines(result, file)
    message(sprintf("LaTeX table written to: %s", file))
    return(invisible(result))
  }

  result
}
