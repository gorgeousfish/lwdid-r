#' @title Warning Registry
#' @description Buffered warning collection for staggered DiD estimation.
#'   Collects warnings during (cohort, period) loops and flushes aggregated
#'   summaries after completion, avoiding warning floods.
#' @name warning-registry
#' @family warning-registry
NULL

# Critical warning categories always shown even in quiet mode
.CRITICAL_CATEGORIES <- c("lwdid_convergence", "lwdid_numerical")

# Valid verbose levels
.VALID_VERBOSE_LEVELS <- c("quiet", "default", "verbose")

#' Create a new warning registry
#'
#' @return An environment-based warning registry object with methods:
#'   `$register()`, `$flush()`, `$get_log()`, `$get_diagnostics()`,
#'   `$count()`, `$clear()`, `$merge()`.
#' @export
#' @family warning-registry
#' @examples
#' reg <- new_warning_registry()
#' reg$register("lwdid_small_sample", "Sample size below 30",
#'              cohort = 2005L, period = 2003L,
#'              context = list(n_treated = 3))
#' reg$count()
new_warning_registry <- function() {
  self <- new.env(parent = emptyenv())
  self$records <- list()
  self$flushed <- FALSE

  self$register <- function(category, message, cohort = NULL, period = NULL,
                            context = list()) {
    # Validate category
    if (!is.character(category) || length(category) != 1L) {
      stop_lwdid(
        sprintf("'category' must be a single character string, got %s of length %d",
                class(category)[1], length(category)),
        class = "lwdid_invalid_parameter",
        param = "category", value = category,
        allowed = "single character string"
      )
    }
    # Validate message
    if (!is.character(message) || length(message) != 1L) {
      stop_lwdid(
        sprintf("'message' must be a single character string, got %s of length %d",
                class(message)[1], length(message)),
        class = "lwdid_invalid_parameter",
        param = "message", value = message,
        allowed = "single character string"
      )
    }
    # Convert NULL context to empty list
    if (is.null(context)) context <- list()

    rec <- list(
      category  = category,
      message   = message,
      cohort    = cohort,
      period    = period,
      context   = context,
      timestamp = as.numeric(Sys.time())
    )
    self$records <- c(self$records, list(rec))
    invisible(NULL)
  }

  self$flush <- function(verbose = "default", total_pairs = NULL) {
    # Step 1: Validate verbose (always, even if already flushed)
    verbose <- tolower(verbose)
    if (!verbose %in% .VALID_VERBOSE_LEVELS) {
      stop_lwdid(
        sprintf("Invalid verbose level '%s'. Allowed: %s",
                verbose, paste(.VALID_VERBOSE_LEVELS, collapse = ", ")),
        class = "lwdid_invalid_parameter",
        param = "verbose", value = verbose,
        allowed = .VALID_VERBOSE_LEVELS
      )
    }
    # Step 2: Check flushed flag
    if (self$flushed) return(invisible(NULL))
    # Step 3: Mark as flushed
    self$flushed <- TRUE
    # Step 4: Check if records empty
    if (length(self$records) == 0L) return(invisible(NULL))
    # Step 5: Dispatch by verbose level
    grouped <- .aggregate_by_category(self$records)

    if (verbose == "verbose") {
      # Output each record's original message individually
      for (rec in self$records) {
        warn_lwdid(rec$message, class = rec$category)
      }
    } else {
      # quiet or default
      cats_to_show <- if (verbose == "quiet") {
        intersect(names(grouped), .CRITICAL_CATEGORIES)
      } else {
        names(grouped)
      }
      for (cat_name in cats_to_show) {
        cat_records <- grouped[[cat_name]]
        summary_msg <- .format_summary(cat_name, cat_records, total_pairs)
        warn_lwdid(summary_msg, class = cat_name)
      }
    }
    invisible(NULL)
  }

  self$get_log <- function() {
    if (length(self$records) == 0L) {
      return(data.frame(
        category  = character(0),
        message   = character(0),
        cohort    = integer(0),
        period    = integer(0),
        timestamp = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      category  = vapply(self$records, `[[`, character(1), "category"),
      message   = vapply(self$records, `[[`, character(1), "message"),
      cohort    = vapply(self$records, function(r) {
        r$cohort %||% NA_integer_
      }, integer(1)),
      period    = vapply(self$records, function(r) {
        r$period %||% NA_integer_
      }, integer(1)),
      timestamp = vapply(self$records, `[[`, numeric(1), "timestamp"),
      stringsAsFactors = FALSE
    )
  }

  self$get_diagnostics <- function() {
    if (length(self$records) == 0L) return(list())
    grouped <- .aggregate_by_category(self$records)
    result <- list()
    for (cat_name in names(grouped)) {
      cat_records <- grouped[[cat_name]]
      # Affected pairs: records where cohort or period is non-NULL
      affected <- Filter(function(r) {
        !is.null(r$cohort) || !is.null(r$period)
      }, cat_records)
      affected_pairs <- lapply(affected, function(r) {
        c(r$cohort %||% NA_integer_, r$period %||% NA_integer_)
      })
      ctx_summary <- .summarize_context(cat_records)
      result <- c(result, list(list(
        category        = cat_name,
        message         = cat_records[[1]]$message,
        count           = length(cat_records),
        affected_pairs  = affected_pairs,
        context_summary = ctx_summary
      )))
    }
    result
  }

  self$count <- function(category = NULL) {
    if (is.null(category)) return(length(self$records))
    as.integer(sum(vapply(self$records, function(r) {
      r$category == category
    }, logical(1))))
  }

  self$clear <- function() {
    self$records <- list()
    self$flushed <- FALSE
    invisible(NULL)
  }

  self$merge <- function(other_registry) {
    # Validate other_registry
    if (!is.environment(other_registry) ||
        is.null(other_registry$records) ||
        is.null(other_registry$flushed)) {
      stop_lwdid(
        "merge() requires a valid warning registry object.",
        class = "lwdid_invalid_parameter",
        param = "other_registry", value = class(other_registry),
        allowed = "warning registry environment"
      )
    }
    self$records <- c(self$records, other_registry$records)
    self$flushed <- FALSE
    invisible(NULL)
  }

  self
}

# ============================================================================
# Internal Helper Functions
# ============================================================================

#' Aggregate records by category preserving insertion order
#' @keywords internal
.aggregate_by_category <- function(records) {
  result <- list()
  seen_order <- character(0)
  for (rec in records) {
    cat_name <- rec$category
    if (!cat_name %in% seen_order) {
      seen_order <- c(seen_order, cat_name)
      result[[cat_name]] <- list()
    }
    result[[cat_name]] <- c(result[[cat_name]], list(rec))
  }
  result
}

#' Format aggregated summary message for a category
#' @keywords internal
.format_summary <- function(cat_name, cat_records, total_pairs = NULL) {
  representative_msg <- cat_records[[1]]$message
  # Count affected pairs (records with non-NULL cohort or period)
  affected <- Filter(function(r) {
    !is.null(r$cohort) || !is.null(r$period)
  }, cat_records)
  n_affected <- length(affected)

  if (n_affected > 0L) {
    if (!is.null(total_pairs) && total_pairs > 0L) {
      sprintf("[%s] %s (%d/%d (cohort, period) pairs)",
              cat_name, representative_msg, n_affected, total_pairs)
    } else {
      sprintf("[%s] %s (%d (cohort, period) pairs)",
              cat_name, representative_msg, n_affected)
    }
  } else {
    sprintf("[%s] %s (%d occurrences)",
            cat_name, representative_msg, length(cat_records))
  }
}

#' Summarize numeric context fields across records
#' @keywords internal
.summarize_context <- function(cat_records) {
  # Collect all context keys
  all_keys <- unique(unlist(lapply(cat_records, function(r) names(r$context))))
  if (length(all_keys) == 0L) return(list())

  # Sort alphabetically
  all_keys <- sort(all_keys)
  result <- list()
  for (key in all_keys) {
    vals <- lapply(cat_records, function(r) r$context[[key]])
    vals <- vals[!vapply(vals, is.null, logical(1))]
    # Only numeric fields
    numeric_vals <- Filter(is.numeric, vals)
    if (length(numeric_vals) > 0L) {
      numeric_vals <- unlist(numeric_vals)
      result[[paste0(key, "_min")]] <- min(numeric_vals, na.rm = TRUE)
      result[[paste0(key, "_max")]] <- max(numeric_vals, na.rm = TRUE)
    }
  }
  result
}
