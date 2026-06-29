#' @title Compare Multiple lwdid Estimation Results
#' @description
#' Compare lwdid estimation results across different specifications
#' (estimators, transformations, control group strategies) in a standardized
#' table format. Based on ATT definitions from Lee & Wooldridge (2026).
#'
#' Different rolling methods estimate the same target parameter tau under
#' different identification assumptions. The mathematical purpose of compare()
#' is to demonstrate ATT robustness across assumptions and methods.
#'
#' @param ... Named lwdid_result objects, or a list containing them
#' @param type Character. Comparison type: "overall" (default, overall ATT) or
#'   "effects" (period-specific effects)
#' @param stats Character vector. Statistics to include. Default:
#'   c("estimate", "std.error", "ci", "p.value")
#' @param stars Logical. Add significance stars? Default TRUE
#' @param digits Integer. Number of decimal places. Default 3
#'
#' @return An \code{lwdid_comparison} S3 object (inherits data.frame), containing:
#'   \itemize{
#'     \item Each model specification as a column
#'     \item Statistics (estimate, std. error, CI, p-value, etc.) as rows
#'     \item Formatted print method
#'   }
#'
#' @seealso \code{\link{tidy.lwdid_result}}, \code{\link{print.lwdid_comparison}}
#' @family lwdid-results
#'
#' @export
#' @examples
#' \donttest{
#'   res1 <- structure(list(
#'     att = 0.5, se_att = 0.1, t_stat = 5.0, pvalue = 0.001,
#'     ci_lower = 0.3, ci_upper = 0.7, nobs = 100L,
#'     n_treated = 50L, n_control = 50L, rolling = "demean",
#'     estimator = "ra", vce_type = "HC3", method = "proc21"
#'   ), class = "lwdid_result")
#'
#'   res2 <- structure(list(
#'     att = 0.6, se_att = 0.15, t_stat = 4.0, pvalue = 0.005,
#'     ci_lower = 0.31, ci_upper = 0.89, nobs = 100L,
#'     n_treated = 50L, n_control = 50L, rolling = "demean",
#'     estimator = "ipw", vce_type = "HC3", method = "proc21"
#'   ), class = "lwdid_result")
#'
#'   # Compare two specifications (overall ATT)
#'   compare(RA = res1, IPW = res2)
#'
#'   # type = "effects" requires results from staggered estimation
#'   # with cohort-level effects (att_by_cohort_time field):
#'   # compare(RA = stag_res1, IPW = stag_res2, type = "effects")
#' }
compare <- function(..., type = c("overall", "effects"),
                    stats = c("estimate", "std.error", "ci", "p.value"),
                    stars = TRUE, digits = 3L) {

  # --- 0. Argument matching ---
  type <- match.arg(type)
  digits <- as.integer(digits)
  stopifnot(is.logical(stars), length(stars) == 1L)
  stopifnot(is.numeric(digits), length(digits) == 1L, digits >= 0L)

  # --- 1. Collect model objects ---
  dots <- list(...)

  # Support passing a single list

  if (length(dots) == 1L && is.list(dots[[1L]]) &&
      !inherits(dots[[1L]], "lwdid_result")) {
    dots <- dots[[1L]]
  }

  n_models <- length(dots)
  if (n_models == 0L) {
    stop("At least one lwdid_result object is required for comparison", call. = FALSE)
  }

  # --- 2. Input validation: ensure all objects are lwdid_result ---
  for (i in seq_len(n_models)) {
    if (!inherits(dots[[i]], "lwdid_result")) {
      stop(sprintf("Argument %d is not a lwdid_result object", i),
           call. = FALSE)
    }
  }

  # --- 3. Auto-naming ---
  model_names <- names(dots)
  if (is.null(model_names) || any(!nzchar(model_names))) {
    # Auto-generate descriptive names for unnamed models
    auto_names <- vapply(dots, function(m) {
      rolling_part <- if (!is.null(m$rolling) && !is.na(m$rolling)) {
        as.character(m$rolling)
      } else {
        "unknown"
      }
      est_part <- if (!is.null(m$estimator) && !is.na(m$estimator)) {
        as.character(m$estimator)
      } else {
        "est"
      }
      paste0(rolling_part, "_", est_part)
    }, character(1))

    if (is.null(model_names)) {
      model_names <- auto_names
    } else {
      # Fill in empty names only
      empty_idx <- !nzchar(model_names)
      model_names[empty_idx] <- auto_names[empty_idx]
    }

    # Handle duplicate names
    if (anyDuplicated(model_names)) {
      model_names <- make.unique(model_names, sep = "_")
    }
  }

  # --- 4. Extract statistics by type ---
  if (type == "overall") {
    result <- .compare_overall(dots, model_names, stats, stars, digits)
  } else {
    result <- .compare_effects(dots, model_names, stats, stars, digits)
  }

  result
}


# =============================================================================
# Internal function: compare overall ATT
# =============================================================================
.compare_overall <- function(models, model_names, stats, stars, digits) {
  n_models <- length(models)

  # Extract core statistics from each model
  att_vals <- vapply(models, function(m) {
    if (!is.null(m$att) && is.finite(m$att)) m$att else NA_real_
  }, numeric(1))

  se_vals <- vapply(models, function(m) {
    if (!is.null(m$se_att) && is.finite(m$se_att)) m$se_att else NA_real_
  }, numeric(1))

  pval_vals <- vapply(models, function(m) {
    if (!is.null(m$pvalue) && is.finite(m$pvalue)) m$pvalue else NA_real_
  }, numeric(1))

  ci_lower_vals <- vapply(models, function(m) {
    if (!is.null(m$ci_lower) && is.finite(m$ci_lower)) m$ci_lower else NA_real_
  }, numeric(1))

  ci_upper_vals <- vapply(models, function(m) {
    if (!is.null(m$ci_upper) && is.finite(m$ci_upper)) m$ci_upper else NA_real_
  }, numeric(1))

  # Generate significance stars
  star_strs <- if (stars) {
    vapply(pval_vals, function(p) {
      if (is.na(p)) return("")
      if (p < 0.01) return("***")
      if (p < 0.05) return("**")
      if (p < 0.10) return("*")
      ""
    }, character(1))
  } else {
    rep("", n_models)
  }

  # Build output table rows
  rows <- list()
  row_labels <- character(0)

  if ("estimate" %in% stats) {
    # ATT estimate row (with stars)
    est_row <- vapply(seq_len(n_models), function(i) {
      if (is.na(att_vals[i])) return("NA")
      paste0(formatC(att_vals[i], format = "f", digits = digits), star_strs[i])
    }, character(1))
    rows <- c(rows, list(est_row))
    row_labels <- c(row_labels, "ATT")
  }

  if ("std.error" %in% stats) {
    # Standard error row (parenthesized)
    se_row <- vapply(se_vals, function(s) {
      if (is.na(s)) return("")
      sprintf("(%s)", formatC(s, format = "f", digits = digits))
    }, character(1))
    rows <- c(rows, list(se_row))
    row_labels <- c(row_labels, "")
  }

  if ("ci" %in% stats) {
    # Confidence interval row
    alpha_val <- models[[1]]$alpha %||% 0.05
    ci_level_pct <- round((1 - alpha_val) * 100)
    ci_row <- vapply(seq_len(n_models), function(i) {
      if (is.na(ci_lower_vals[i]) || is.na(ci_upper_vals[i])) return("")
      sprintf("[%s,%s]",
              formatC(ci_lower_vals[i], format = "f", digits = digits),
              formatC(ci_upper_vals[i], format = "f", digits = digits))
    }, character(1))
    rows <- c(rows, list(ci_row))
    row_labels <- c(row_labels, sprintf("%d%% CI", ci_level_pct))
  }

  if ("p.value" %in% stats) {
    # p-value row
    pval_row <- vapply(pval_vals, function(p) {
      if (is.na(p)) return("")
      if (p < 0.001) return("<0.001")
      formatC(p, format = "f", digits = digits)
    }, character(1))
    rows <- c(rows, list(pval_row))
    row_labels <- c(row_labels, "p-value")
  }

  # Model info rows
  info_rows <- list()
  info_labels <- character(0)

  # N (number of observations)
  n_row <- vapply(models, function(m) {
    n <- m$nobs
    if (!is.null(n) && is.finite(n)) as.character(as.integer(n)) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(n_row))
  info_labels <- c(info_labels, "N")

  # N_treated
  nt_row <- vapply(models, function(m) {
    n <- m$n_treated
    if (!is.null(n) && is.finite(n)) as.character(as.integer(n)) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(nt_row))
  info_labels <- c(info_labels, "N_treated")

  # N_control
  nc_row <- vapply(models, function(m) {
    n <- m$n_control
    if (!is.null(n) && is.finite(n)) as.character(as.integer(n)) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(nc_row))
  info_labels <- c(info_labels, "N_control")

  # Transformation (rolling method)
  trans_row <- vapply(models, function(m) {
    if (!is.null(m$rolling) && !is.na(m$rolling)) as.character(m$rolling) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(trans_row))
  info_labels <- c(info_labels, "Transformation")

  # Estimator
  est_info_row <- vapply(models, function(m) {
    if (!is.null(m$estimator) && !is.na(m$estimator)) as.character(m$estimator) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(est_info_row))
  info_labels <- c(info_labels, "Estimator")

  # VCE
  vce_row <- vapply(models, function(m) {
    if (!is.null(m$vce_type) && !is.na(m$vce_type)) as.character(m$vce_type) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(vce_row))
  info_labels <- c(info_labels, "VCE")

  # Control group (strategy; only shown if not all identical)
  cg_vals <- vapply(models, function(m) {
    if (!is.null(m$control_group) && !is.na(m$control_group)) {
      as.character(m$control_group)
    } else {
      "NA"
    }
  }, character(1))
  if (length(unique(cg_vals)) > 1L) {
    info_rows <- c(info_rows, list(cg_vals))
    info_labels <- c(info_labels, "Control group")
  }

  # Assemble result data.frame
  all_rows <- c(rows, info_rows)
  all_labels <- c(row_labels, info_labels)

  # Build data.frame
  result_df <- as.data.frame(
    matrix(unlist(all_rows), nrow = length(all_rows), byrow = TRUE),
    stringsAsFactors = FALSE
  )
  names(result_df) <- model_names
  result_df <- cbind(data.frame(stat = all_labels, stringsAsFactors = FALSE),
                     result_df)
  rownames(result_df) <- NULL

  # Attach metadata
  attr(result_df, "n_models") <- n_models
  attr(result_df, "model_names") <- model_names
  attr(result_df, "type") <- "overall"
  attr(result_df, "stars") <- stars
  attr(result_df, "digits") <- digits
  attr(result_df, "n_stat_rows") <- length(rows)

  class(result_df) <- c("lwdid_comparison", "data.frame")
  result_df
}


# =============================================================================
# Internal function: compare effects (period-specific)
# =============================================================================
.compare_effects <- function(models, model_names, stats, stars, digits) {
  n_models <- length(models)

  # Try to extract period-specific effects from effects field

  effects_list <- lapply(models, function(m) {
    eff <- m$effects
    if (is.null(eff) || !is.data.frame(eff) || nrow(eff) == 0L) {
      eff <- m$att_by_cohort_time
    }
    if (is.null(eff) || !is.data.frame(eff) || nrow(eff) == 0L) {
      eff <- m$cohort_effects
    }
    if (is.null(eff) || !is.data.frame(eff) || nrow(eff) == 0L) {
      eff <- m$event_time_effects
    }
    eff
  })

  # Check if at least one model has period-specific effects
  has_effects <- vapply(effects_list, function(e) {
    !is.null(e) && is.data.frame(e) && nrow(e) > 0L
  }, logical(1))

  if (!any(has_effects)) {
    message("No models contain period-specific effects; falling back to type='overall'.")
    return(.compare_overall(models, model_names, stats, stars, digits))
  }

  # Align effects across models by (cohort, period) or (event_time) keys
  # to avoid silently comparing different rows positionally
  ref_idx <- which(has_effects)[1L]
  ref_eff <- effects_list[[ref_idx]]

  # Determine identifier columns for row alignment
  id_cols <- intersect(names(ref_eff), c("cohort", "period", "event_time", "g", "r", "time"))

  if (length(id_cols) >= 1L) {
    # Build a unified key set from ALL models (union of all cohort-period combos)
    all_keys <- do.call(rbind, lapply(which(has_effects), function(i) {
      effects_list[[i]][, id_cols, drop = FALSE]
    }))
    all_keys <- unique(all_keys)
    all_keys <- all_keys[do.call(order, as.list(all_keys)), , drop = FALSE]
    rownames(all_keys) <- NULL

    # Reindex each model's effects to match the unified key order
    effects_list <- lapply(effects_list, function(eff) {
      if (is.null(eff) || !is.data.frame(eff) || nrow(eff) == 0L) return(NULL)
      merged <- merge(all_keys, eff, by = id_cols, all.x = TRUE, sort = FALSE)
      # Restore unified order
      merged <- merged[do.call(order, as.list(merged[, id_cols, drop = FALSE])), , drop = FALSE]
      rownames(merged) <- NULL
      merged
    })

    period_labels <- apply(all_keys, 1, function(row) paste(row, collapse = ","))
  } else {
    # No id columns: fall back to positional (row number) comparison
    period_labels <- paste0("Period_", seq_len(nrow(ref_eff)))
  }

  n_periods <- length(period_labels)

  # Extract ATT and SE for each model/period
  result_rows <- list()
  result_labels <- character(0)

  for (p in seq_len(n_periods)) {
    # ATT 行
    att_row <- vapply(seq_len(n_models), function(i) {
      eff <- effects_list[[i]]
      if (is.null(eff) || nrow(eff) < p) return("NA")
      att_val <- eff$att[p]
      pval <- if ("pvalue" %in% names(eff)) eff$pvalue[p] else NA_real_
      star_str <- ""
      if (stars && !is.na(pval)) {
        if (pval < 0.01) star_str <- "***"
        else if (pval < 0.05) star_str <- "**"
        else if (pval < 0.10) star_str <- "*"
      }
      if (is.na(att_val)) return("NA")
      paste0(formatC(att_val, format = "f", digits = digits), star_str)
    }, character(1))
    result_rows <- c(result_rows, list(att_row))
    result_labels <- c(result_labels, period_labels[p])

    # SE 行（如果需要）
    if ("std.error" %in% stats) {
      se_row <- vapply(seq_len(n_models), function(i) {
        eff <- effects_list[[i]]
        if (is.null(eff) || nrow(eff) < p) return("")
        se_val <- if ("se" %in% names(eff)) eff$se[p] else NA_real_
        if (is.na(se_val)) return("")
        sprintf("(%s)", formatC(se_val, format = "f", digits = digits))
      }, character(1))
      result_rows <- c(result_rows, list(se_row))
      result_labels <- c(result_labels, "")
    }
  }

  # Assemble data.frame
  result_df <- as.data.frame(
    matrix(unlist(result_rows), nrow = length(result_rows), byrow = TRUE),
    stringsAsFactors = FALSE
  )
  names(result_df) <- model_names
  result_df <- cbind(data.frame(stat = result_labels, stringsAsFactors = FALSE),
                     result_df)
  rownames(result_df) <- NULL

  # Attach metadata
  attr(result_df, "n_models") <- n_models
  attr(result_df, "model_names") <- model_names
  attr(result_df, "type") <- "effects"
  attr(result_df, "stars") <- stars
  attr(result_df, "digits") <- digits
  attr(result_df, "n_periods") <- n_periods

  class(result_df) <- c("lwdid_comparison", "data.frame")
  result_df
}


# =============================================================================
# S3 method: formatted printing
# =============================================================================

#' Print method for lwdid_comparison objects
#'
#' @description
#' Formats and prints a comparison table showing side-by-side results
#' from multiple lwdid specifications.
#'
#' @param x An lwdid_comparison object returned by \code{\link{compare}}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#'
#' @seealso \code{\link{compare}}
#' @family lwdid-results
#' @export
print.lwdid_comparison <- function(x, ...) {
  type <- attr(x, "type") %||% "overall"
  model_names <- attr(x, "model_names") %||% names(x)[-1L]
  stars_used <- attr(x, "stars") %||% TRUE
  n_stat_rows <- attr(x, "n_stat_rows") %||% 0L

  # Calculate column widths
  n_models <- length(model_names)
  col_widths <- vapply(seq_len(n_models), function(i) {
    col_data <- x[[i + 1L]]  # skip stat column
    max(nchar(model_names[i]), max(nchar(as.character(col_data)), na.rm = TRUE))
  }, integer(1))
  # Ensure minimum width

  col_widths <- pmax(col_widths, 10L)

  label_width <- max(nchar(x$stat), 16L)

  # Total table width
  total_width <- label_width + 2L + sum(col_widths + 2L)

  # Title
  cat("\nlwdid Comparison Table\n")
  cat(strrep("=", total_width), "\n")

  # Header row
  header <- sprintf("%-*s", label_width, "")
  for (i in seq_len(n_models)) {
    header <- paste0(header, "  ", sprintf("%-*s", col_widths[i], model_names[i]))
  }
  cat(header, "\n")
  cat(strrep("-", total_width), "\n")

  # Statistics rows
  n_rows <- nrow(x)
  for (r in seq_len(n_rows)) {
    row_label <- x$stat[r]

    # Add separator between statistic rows and info rows
    if (n_stat_rows > 0L && r == n_stat_rows + 1L) {
      cat(strrep("-", total_width), "\n")
    }

    line <- sprintf("%-*s", label_width, row_label)
    for (i in seq_len(n_models)) {
      cell <- as.character(x[[i + 1L]][r])
      line <- paste0(line, "  ", sprintf("%-*s", col_widths[i], cell))
    }
    cat(line, "\n")
  }

  # Footer
  cat(strrep("=", total_width), "\n")
  if (stars_used) {
    cat("*** p<0.01, ** p<0.05, * p<0.1\n")
  }
  cat("\n")

  invisible(x)
}
