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
#' @export
to_csv <- function(x, ...) UseMethod("to_csv")

#' Convert an object to a dictionary (named list)
#'
#' @param x An object.
#' @param ... Additional arguments passed to methods.
#' @return A named list of parameters and statistics.
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
  if (!is.null(caption)) lines <- c(lines, sprintf("\\caption{%s}", caption))
  if (!is.null(label)) lines <- c(lines, sprintf("\\label{%s}", label))
  lines <- c(lines, "\\begin{threeparttable}")
  lines <- c(lines, "\\begin{tabular}{lc}")
  lines <- c(lines, rule_top)

  # --- 3. Main body rows ---
  lines <- c(lines, sprintf("ATT & %s%s \\\\",
                              sprintf("%.*f", digits, x$att), star_str))

  if (se_in_parentheses) {
    lines <- c(lines, sprintf(" & (%s) \\\\",
                                sprintf("%.*f", digits, x$se_att)))
  } else {
    lines <- c(lines, sprintf("SE & %s \\\\",
                                sprintf("%.*f", digits, x$se_att)))
  }

  lines <- c(lines, rule_mid)

  if (include_ci) {
    lines <- c(lines, sprintf("CI [%d%%] & [%s, %s] \\\\",
                                ci_level_pct,
                                sprintf("%.*f", digits, x$ci_lower),
                                sprintf("%.*f", digits, x$ci_upper)))
  }

  # p-value row
  if (is.na(x$pvalue)) {
    lines <- c(lines, "p-value & -- \\\\")
  } else if (x$pvalue < 0.001) {
    lines <- c(lines, "p-value & $<$0.001 \\\\")
  } else {
    lines <- c(lines, sprintf("p-value & %s \\\\",
                                sprintf("%.*f", digits, x$pvalue)))
  }

  lines <- c(lines, rule_mid)

  # Sample info rows
  lines <- c(lines, sprintf("N & %d \\\\", x$nobs))
  lines <- c(lines, sprintf("$N_{treated}$ & %d \\\\", x$n_treated))
  lines <- c(lines, sprintf("$N_{control}$ & %d \\\\", x$n_control))
  lines <- c(lines, sprintf("Estimator & %s \\\\", .escape_latex(x$estimator)))
  lines <- c(lines, sprintf("VCE & %s \\\\", .escape_latex(x$vce_type)))
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
                                sprintf("%.*f", digits, x$ri_pvalue)))
    if (!is.null(x$ri_seed)) {
      lines <- c(lines, sprintf("RI seed & %d \\\\", as.integer(x$ri_seed)))
    }
    if (!is.null(x$rireps)) {
      lines <- c(lines, sprintf("RI reps & %d \\\\", as.integer(x$rireps)))
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
    lines <- c(lines, "\\vspace{1em}")
    lines <- c(lines, "\\begin{threeparttable}")
    lines <- c(lines, "\\begin{tabular}{lccccc}")
    lines <- c(lines, rule_top)
    lines <- c(lines, "Period & ATT & SE & t-stat & p-value & CI \\\\")
    lines <- c(lines, rule_mid)

    for (i in seq_len(nrow(x$att_by_period))) {
      row <- x$att_by_period[i, ]

      # Period identifier (period or event_time)
      period_id <- if ("period" %in% names(row)) row$period
                   else if ("event_time" %in% names(row)) row$event_time
                   else i

      # ATT and SE
      att_val <- sprintf("%.*f", digits, row$att)
      se_val <- sprintf("%.*f", digits, row$se)

      # Period stars
      period_star <- ""
      if (stars && "pvalue" %in% names(row) && !is.na(row$pvalue)) {
        if (row$pvalue < 0.01) period_star <- "***"
        else if (row$pvalue < 0.05) period_star <- "**"
        else if (row$pvalue < 0.10) period_star <- "*"
      }

      # t-statistic
      t_val <- if ("t_stat" %in% names(row)) sprintf("%.*f", digits, row$t_stat)
               else if ("t" %in% names(row)) sprintf("%.*f", digits, row$t)
               else sprintf("%.*f", digits, row$att / row$se)

      # p-value
      p_str <- if (!"pvalue" %in% names(row) || is.na(row$pvalue)) "--"
               else if (row$pvalue < 0.001) "$<$0.001"
               else sprintf("%.*f", digits, row$pvalue)

      # CI (handle column name variants)
      ci_lo <- if ("ci_lower" %in% names(row)) row$ci_lower
               else if ("conf.low" %in% names(row)) row[["conf.low"]]
               else NA
      ci_hi <- if ("ci_upper" %in% names(row)) row$ci_upper
               else if ("conf.high" %in% names(row)) row[["conf.high"]]
               else NA
      ci_str <- if (is.na(ci_lo) || is.na(ci_hi)) "--"
                else sprintf("[%s, %s]",
                             sprintf("%.*f", digits, ci_lo),
                             sprintf("%.*f", digits, ci_hi))

      lines <- c(lines, sprintf("%s & %s%s & %s & %s & %s & %s \\\\",
                                  period_id, att_val, period_star,
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
    lines <- c(lines, "\\vspace{1em}")
    lines <- c(lines, "\\begin{threeparttable}")
    lines <- c(lines, "\\begin{tabular}{lcccc}")
    lines <- c(lines, rule_top)
    lines <- c(lines, "Cohort & ATT & SE & N & Weight \\\\")
    lines <- c(lines, rule_mid)

    for (i in seq_len(nrow(x$att_by_cohort))) {
      row <- x$att_by_cohort[i, ]
      cohort_id <- if ("cohort" %in% names(row)) row$cohort else i
      cohort_att <- sprintf("%.*f", digits, row$att)
      cohort_se <- sprintf("%.*f", digits, row$se)
      cohort_n <- if ("n" %in% names(row)) sprintf("%d", row$n)
                  else if ("nobs" %in% names(row)) sprintf("%d", row$nobs)
                  else "--"

      # Weight: prefer att_by_cohort$weight, fallback to cohort_weights
      cohort_w <- if ("weight" %in% names(row)) {
        sprintf("%.*f", digits, row$weight)
      } else if (!is.null(x$cohort_weights)) {
        cid <- as.character(cohort_id)
        if (cid %in% names(x$cohort_weights)) {
          sprintf("%.*f", digits, x$cohort_weights[[cid]])
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
                                  cohort_id, cohort_att, cohort_star,
                                  cohort_se, cohort_n, cohort_w))
    }

    lines <- c(lines, rule_mid)
    lines <- c(lines, sprintf("Overall (weighted) & %s%s & %s & %d & -- \\\\",
                                sprintf("%.*f", digits, x$att), star_str,
                                sprintf("%.*f", digits, x$se_att),
                                x$nobs))
    lines <- c(lines, rule_bot)
    lines <- c(lines, "\\end{tabular}")
    lines <- c(lines, "\\end{threeparttable}")
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

#' @export
to_csv.lwdid_result <- function(x, file, what = "summary", ...) {
  # --- 0. Parameter validation ---
  stopifnot(is.character(file), length(file) == 1L)
  what <- match.arg(what, c("summary", "by_period", "by_cohort", "all"))

  # --- 1. Build data.frame based on what ---
  df <- switch(what,
    "summary" = {
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
        rolling = as.character(x$rolling %||% NA_character_),
        stringsAsFactors = FALSE
      )
      if (!is.null(x$ri_pvalue)) {
        base_df$ri_pvalue <- x$ri_pvalue
        base_df$ri_seed <- if (!is.null(x$ri_seed)) x$ri_seed else NA
        base_df$rireps <- if (!is.null(x$rireps)) x$rireps else NA
      }
      base_df
    },
    "by_period" = {
      if (is.null(x$att_by_period) ||
          (is.data.frame(x$att_by_period) && nrow(x$att_by_period) == 0L)) {
        stop("No period-specific results available")
      }
      x$att_by_period
    },
    "by_cohort" = {
      if (is.null(x$att_by_cohort) ||
          (is.data.frame(x$att_by_cohort) && nrow(x$att_by_cohort) == 0L)) {
        stop("No cohort-specific results (Staggered mode only).",
             call. = FALSE)
      }
      x$att_by_cohort
    },
    "all" = {
      if (!is.null(x$att_by_cohort_time)) x$att_by_cohort_time
      else if (!is.null(x$att_by_period)) x$att_by_period
      else data.frame(att = x$att, se = x$se_att)
    }
  )

  # --- 2. Write CSV ---
  utils::write.csv(df, file = file, row.names = FALSE, ...)
  message(sprintf("CSV written to: %s (%d rows)", file, nrow(df)))
  invisible(df)
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

  rule_top <- if (booktabs) "\\toprule" else "\\hline"
  rule_mid <- if (booktabs) "\\midrule" else "\\hline"
  rule_bot <- if (booktabs) "\\bottomrule" else "\\hline"

  col_spec <- paste0("l", paste(rep("c", n_models), collapse = ""))

  lines <- character(0)

  # --- 3. Table header ---
  lines <- c(lines, "\\begin{table}[htbp]")
  lines <- c(lines, "\\centering")
  if (!is.null(caption)) lines <- c(lines, sprintf("\\caption{%s}", caption))
  if (!is.null(label)) lines <- c(lines, sprintf("\\label{%s}", label))
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
    sprintf("%s%s", sprintf("%.*f", digits, models[[i]]$att), star_strs[i])
  }, character(1))
  lines <- c(lines, sprintf("ATT & %s \\\\", paste(att_cells, collapse = " & ")))

  # --- 6. SE row (always parentheses) ---
  se_cells <- vapply(models, function(m) {
    sprintf("(%s)", sprintf("%.*f", digits, m$se_att))
  }, character(1))
  lines <- c(lines, sprintf(" & %s \\\\", paste(se_cells, collapse = " & ")))

  # --- 7. CI row (optional) ---
  if (include_ci) {
    ci_cells <- vapply(models, function(m) {
      ci_level_pct <- round((1 - m$alpha) * 100)
      sprintf("[%s, %s]",
              sprintf("%.*f", digits, m$ci_lower),
              sprintf("%.*f", digits, m$ci_upper))
    }, character(1))
    ci_level_pct <- round((1 - models[[1]]$alpha) * 100)
    lines <- c(lines, sprintf("CI [%d%%] & %s \\\\",
                                ci_level_pct,
                                paste(ci_cells, collapse = " & ")))
  }

  # --- 8. p-value row ---
  pval_cells <- vapply(models, function(m) {
    if (is.na(m$pvalue)) return("--")
    if (m$pvalue < 0.001) return("$<$0.001")
    sprintf("%.*f", digits, m$pvalue)
  }, character(1))
  lines <- c(lines, sprintf("p-value & %s \\\\",
                              paste(pval_cells, collapse = " & ")))
  lines <- c(lines, rule_mid)

  # --- 9. Sample info rows ---
  n_cells <- vapply(models, function(m) sprintf("%d", m$nobs), character(1))
  lines <- c(lines, sprintf("N & %s \\\\", paste(n_cells, collapse = " & ")))

  est_cells <- vapply(models, function(m) .escape_latex(m$estimator), character(1))
  lines <- c(lines, sprintf("Estimator & %s \\\\",
                              paste(est_cells, collapse = " & ")))

  vce_cells <- vapply(models, function(m) .escape_latex(m$vce_type), character(1))
  lines <- c(lines, sprintf("VCE & %s \\\\",
                              paste(vce_cells, collapse = " & ")))

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
