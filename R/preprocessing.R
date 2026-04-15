# ============================================================================
# preprocessing.R — Data Preprocessing Module for lwdid
#
# Implements aggregate_to_panel() for collapsing repeated cross-sectional data
# to panel format, following Lee & Wooldridge (2026) Section 3.
#
# Formula: Y_bar_st = sum_{i in (s,t)} w_ist * Y_ist, where sum w_ist = 1
#
# Structure:
#   - Story 11.1: Input validation functions
#   - Story 11.2: Weight normalization & weighted computation core
#   - Story 11.5: CellStatistics constructor
#   - Story 11.4: AggregationResult S3 class & output methods
#   - Story 11.3: Main aggregate_to_panel() function & helpers
# ============================================================================

# ============================================================================
# Constants
# ============================================================================

WEIGHT_SUM_TOLERANCE <- 1e-9

# ============================================================================
# Story 11.1: Input Validation Functions
# ============================================================================

#' Validate aggregation input data and column existence
#' @keywords internal
.validate_aggregation_inputs <- function(data, unit_var, time_var,
                                          outcome_var, weight_var = NULL,
                                          controls = NULL) {
  # Type check
  if (!is.data.frame(data)) {
    stop_lwdid(
      sprintf("Input data must be a data.frame or data.table. Got: %s",
              class(data)[1]),
      class = "lwdid_invalid_parameter",
      param = "data", value = class(data)[1],
      allowed = c("data.frame", "data.table")
    )
  }

  # Empty data check
  if (nrow(data) == 0L) {
    stop_lwdid("Input data is empty",
               class = "lwdid_invalid_parameter",
               param = "data", value = "0 rows")
  }

  # Build required columns list
  required_cols <- c(unit_var, outcome_var)
  required_cols <- c(required_cols, time_var)
  if (!is.null(weight_var)) required_cols <- c(required_cols, weight_var)
  if (!is.null(controls)) required_cols <- c(required_cols, controls)

  # Missing columns check
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop_lwdid(
      sprintf("Required column(s) not found in data: %s. Available columns: %s",
              paste(missing_cols, collapse = ", "),
              paste(names(data), collapse = ", ")),
      class = "lwdid_missing_column",
      column = missing_cols, available = names(data)
    )
  }

  # Outcome variable numeric type check
  if (!is.numeric(data[[outcome_var]])) {
    stop_lwdid(
      sprintf("Outcome variable '%s' must be numeric type. Found: '%s'",
              outcome_var, class(data[[outcome_var]])[1]),
      class = "lwdid_invalid_parameter",
      param = "outcome_var", value = class(data[[outcome_var]])[1]
    )
  }

  invisible(NULL)
}


#' Validate weights: non-negative, handle missing
#' @keywords internal
.validate_weights <- function(data, weight_var) {
  weights <- data[[weight_var]]

  # Negative weight check (before NA handling)
  negative_mask <- weights < 0 & !is.na(weights)
  if (any(negative_mask)) {
    n_negative <- sum(negative_mask)
    examples <- head(weights[negative_mask], 5L)
    stop_lwdid(
      sprintf("Weights must be non-negative. Found %d negative values. Examples: %s",
              n_negative, paste(round(examples, 6), collapse = ", ")),
      class = "lwdid_invalid_parameter",
      param = weight_var, value = n_negative
    )
  }

  # Missing weight handling
  missing_mask <- is.na(weights)
  n_missing <- sum(missing_mask)

  if (n_missing > 0L) {
    warn_lwdid(
      sprintf("Excluded %d observations with missing weights.", n_missing),
      class = "lwdid_data",
      detail = "missing_weights", action_taken = "rows_removed"
    )
    data_clean <- data[!missing_mask]
  } else {
    data_clean <- data.table::copy(data)
  }

  list(data = data_clean, n_missing = n_missing)
}


#' Validate treatment consistency within cells
#' @keywords internal
.validate_treatment_consistency <- function(data, unit_var, time_var,
                                             treatment_var) {
  group_cols <- c(unit_var, time_var)

  within_cell_sd <- data[, .(sd_treat = sd(.SD[[1L]], na.rm = TRUE)),
                          by = group_cols, .SDcols = treatment_var]

  varying <- within_cell_sd[!is.na(sd_treat) & sd_treat > 1e-10]

  if (nrow(varying) > 0L) {
    examples <- head(varying, 3L)
    example_strs <- apply(
      examples[, ..group_cols], 1L,
      function(row) paste(row, collapse = ", ")
    )
    stop_lwdid(
      sprintf(
        paste0("Treatment status varies within %d cell(s). ",
               "Treatment must be constant within each (unit, period) cell. ",
               "Examples of cells with varying treatment: %s"),
        nrow(varying), paste(example_strs, collapse = "; ")
      ),
      class = c("lwdid_invalid_aggregation", "lwdid_aggregation_error"),
      aggregate_type = "treatment_consistency",
      violation = "within_cell_variation",
      detail = sprintf("%d cells with varying treatment", nrow(varying))
    )
  }

  invisible(NULL)
}


#' Validate gvar consistency within units
#' @keywords internal
.validate_gvar_consistency <- function(data, unit_var, gvar) {
  gvar_col <- gvar
  within_unit_sd <- data[, .(sd_gvar = sd(.SD[[1L]], na.rm = TRUE)),
                          by = unit_var, .SDcols = gvar_col]

  varying <- within_unit_sd[!is.na(sd_gvar) & sd_gvar > 1e-10]

  if (nrow(varying) > 0L) {
    examples <- head(varying[[unit_var]], 3L)
    stop_lwdid(
      sprintf(
        paste0("Treatment timing (gvar) varies within %d unit(s). ",
               "gvar must be constant within each unit across all periods. ",
               "Examples of units with varying gvar: %s"),
        nrow(varying), paste(examples, collapse = ", ")
      ),
      class = c("lwdid_invalid_aggregation", "lwdid_aggregation_error"),
      aggregate_type = "gvar_consistency",
      violation = "within_unit_variation",
      detail = sprintf("%d units with varying gvar", nrow(varying))
    )
  }

  invisible(NULL)
}


# ============================================================================
# Story 11.2: Weight Normalization & Weighted Computation Core
# ============================================================================

#' Kahan compensated summation (internal)
#' @param x numeric vector
#' @return compensated sum
#' @keywords internal
.kahan_sum <- function(x) {
  if (length(x) == 0L) return(0.0)
  sum_val <- 0.0
  compensation <- 0.0
  for (i in seq_along(x)) {
    y <- x[i] - compensation
    t <- sum_val + y
    compensation <- (t - sum_val) - y
    sum_val <- t
  }
  sum_val
}


#' Normalize weights to sum to 1
#' @keywords internal
.normalize_weights <- function(weights) {
  if (length(weights) == 0L) return(numeric(0))

  weight_sum <- .kahan_sum(weights)

  if (abs(weight_sum) < 1e-15) {
    n <- length(weights)
    return(rep(1.0 / n, n))
  }

  weights / weight_sum
}


#' Compute weighted mean and optionally variance for a single cell
#' @keywords internal
.compute_cell_weighted_mean <- function(outcomes, raw_weights = NULL,
                                         compute_variance = FALSE) {
  n_total <- length(outcomes)
  valid_idx <- which(!is.na(outcomes))
  outcomes <- outcomes[valid_idx]
  n_obs <- length(outcomes)

  if (n_obs == 0L) {
    return(list(mean = NaN, variance = NULL, ess = NULL, n_obs = 0L))
  }

  # Normalize weights
  if (!is.null(raw_weights)) {
    if (length(raw_weights) != n_total) {
      stop_lwdid(
        sprintf("outcomes length (%d) != raw_weights length (%d)",
                n_total, length(raw_weights)),
        class = "lwdid_invalid_parameter"
      )
    }
    raw_weights <- raw_weights[valid_idx]
    weights <- .normalize_weights(raw_weights)
  } else {
    weights <- rep(1.0 / n_obs, n_obs)
  }

  # Weighted mean
  mean_val <- .kahan_sum(weights * outcomes)

  # Weighted variance (optional)
  variance <- NULL
  if (compute_variance && n_obs > 1L) {
    deviations_sq <- (outcomes - mean_val)^2
    variance <- .kahan_sum(weights * deviations_sq)
  }

  # ESS (survey weights only)
  ess <- NULL
  if (!is.null(raw_weights)) {
    ess <- .compute_effective_sample_size(raw_weights)
  }

  list(mean = mean_val, variance = variance,
       ess = ess, n_obs = n_obs)
}


#' Compute effective sample size
#' @keywords internal
.compute_effective_sample_size <- function(weights) {
  if (length(weights) == 0L) return(0.0)

  weight_sum <- .kahan_sum(weights)
  weight_sq_sum <- .kahan_sum(weights^2)

  if (abs(weight_sq_sum) < 1e-15) {
    return(0.0)
  }

  (weight_sum^2) / weight_sq_sum
}


#' Compute weighted mean for a control variable
#' @keywords internal
.compute_control_weighted_mean <- function(values, raw_weights = NULL) {
  n_total <- length(values)
  valid_idx <- which(!is.na(values))
  values <- values[valid_idx]
  n_obs <- length(values)

  if (n_obs == 0L) return(NaN)

  if (!is.null(raw_weights)) {
    if (length(raw_weights) != n_total) {
      stop_lwdid(
        sprintf("values length (%d) != raw_weights length (%d)",
                n_total, length(raw_weights)),
        class = "lwdid_invalid_parameter"
      )
    }
    raw_weights <- raw_weights[valid_idx]
    weights <- .normalize_weights(raw_weights)
  } else {
    weights <- rep(1.0 / n_obs, n_obs)
  }

  .kahan_sum(weights * values)
}


# ============================================================================
# Story 11.5: CellStatistics Constructor
# ============================================================================

#' Create a CellStatistics named list
#'
#' Converts NULL fields to NA_real_ for rbindlist() compatibility.
#'
#' @param unit Unit identifier value
#' @param period Period identifier value
#' @param n_obs Number of observations in the cell
#' @param outcome_mean Weighted mean of the outcome
#' @param outcome_variance Weighted variance (NULL if not computed)
#' @param effective_sample_size ESS (NULL for equal weights)
#' @param weight_type "equal" or "survey"
#' @return Named list with CellStatistics fields
#' @keywords internal
new_cell_statistics <- function(unit, period, n_obs, outcome_mean,
                                 outcome_variance = NULL,
                                 effective_sample_size = NULL,
                                 weight_type = "equal") {
  list(
    unit = unit,
    period = if (is.list(period)) paste(period, collapse = ", ") else period,
    n_obs = as.integer(n_obs),
    outcome_mean = as.double(outcome_mean),
    outcome_variance = if (is.null(outcome_variance)) NA_real_
                       else as.double(outcome_variance),
    effective_sample_size = if (is.null(effective_sample_size)) NA_real_
                            else as.double(effective_sample_size),
    weight_type = as.character(weight_type)
  )
}


# ============================================================================
# Story 11.4: AggregationResult S3 Class & Output Methods
# ============================================================================

#' Create an lwdid_aggregation_result S3 object
#'
#' @param panel_data data.table, aggregated panel data
#' @param n_original_obs integer, total observations in original data
#' @param n_cells integer, number of (unit, period) cells
#' @param n_units integer, number of unique units
#' @param n_periods integer, number of unique periods
#' @param cell_stats data.table, cell-level statistics
#' @param min_cell_size_stat integer, minimum cell size
#' @param max_cell_size integer, maximum cell size
#' @param mean_cell_size numeric, mean cell size
#' @param median_cell_size numeric, median cell size
#' @param unit_var character, unit variable name
#' @param time_var character, time variable name(s)
#' @param outcome_var character, outcome variable name
#' @param weight_var character or NULL, weight variable name
#' @param frequency character, aggregation frequency
#' @param n_excluded_cells integer, number of excluded cells
#' @param excluded_cells_info list, excluded cell details
#' @return lwdid_aggregation_result S3 object
#' @export
new_lwdid_aggregation_result <- function(
  panel_data, n_original_obs, n_cells, n_units,
  n_periods, cell_stats, min_cell_size_stat,
  max_cell_size, mean_cell_size, median_cell_size,
  unit_var, time_var, outcome_var, weight_var,
  frequency, n_excluded_cells = 0L,
  excluded_cells_info = list()
) {
  if (!data.table::is.data.table(panel_data)) {
    stop("panel_data must be a data.table", call. = FALSE)
  }
  if (!data.table::is.data.table(cell_stats)) {
    stop("cell_stats must be a data.table", call. = FALSE)
  }

  obj <- list(
    panel_data = panel_data,
    n_original_obs = as.integer(n_original_obs),
    n_cells = as.integer(n_cells),
    n_units = as.integer(n_units),
    n_periods = as.integer(n_periods),
    cell_stats = cell_stats,
    min_cell_size_stat = as.integer(min_cell_size_stat),
    max_cell_size = as.integer(max_cell_size),
    mean_cell_size = as.double(mean_cell_size),
    median_cell_size = as.double(median_cell_size),
    unit_var = unit_var,
    time_var = time_var,
    outcome_var = outcome_var,
    weight_var = weight_var,
    frequency = frequency,
    n_excluded_cells = as.integer(n_excluded_cells),
    excluded_cells_info = excluded_cells_info
  )
  class(obj) <- "lwdid_aggregation_result"
  obj
}


#' @export
print.lwdid_aggregation_result <- function(x, ...) {
  cat("Aggregation Summary\n")
  cat("===================\n")
  cat(sprintf("Original observations: %s\n",
              format(x$n_original_obs, big.mark = ",")))
  cat(sprintf("Output cells: %s\n",
              format(x$n_cells, big.mark = ",")))
  cat(sprintf("Units: %s\n",
              format(x$n_units, big.mark = ",")))
  cat(sprintf("Periods: %s\n",
              format(x$n_periods, big.mark = ",")))
  cat(sprintf("Frequency: %s\n", x$frequency))
  wt_str <- if (is.null(x$weight_var)) {
    "None (equal weights)"
  } else {
    x$weight_var
  }
  cat(sprintf("Weight: %s\n", wt_str))
  if (x$n_excluded_cells > 0L) {
    cat(sprintf("Excluded cells: %d\n", x$n_excluded_cells))
  }
  invisible(x)
}


#' @export
summary.lwdid_aggregation_result <- function(object, ...) {
  cat("Aggregation Summary\n")
  cat("===================\n")
  cat(sprintf("Original observations: %s\n",
              format(object$n_original_obs, big.mark = ",")))
  cat(sprintf("Output cells: %s\n",
              format(object$n_cells, big.mark = ",")))
  cat(sprintf("Units: %s\n",
              format(object$n_units, big.mark = ",")))
  cat(sprintf("Periods: %s\n",
              format(object$n_periods, big.mark = ",")))
  cat("\nCell Size Statistics\n")
  cat("--------------------\n")
  cat(sprintf("Minimum: %s\n",
              format(object$min_cell_size_stat, big.mark = ",")))
  cat(sprintf("Maximum: %s\n",
              format(object$max_cell_size, big.mark = ",")))
  cat(sprintf("Mean: %.2f\n", object$mean_cell_size))
  cat(sprintf("Median: %.2f\n", object$median_cell_size))
  cat("\nConfiguration\n")
  cat("-------------\n")
  cat(sprintf("Unit variable: %s\n", object$unit_var))
  tv_str <- paste(object$time_var, collapse = ", ")
  cat(sprintf("Time variable: %s\n", tv_str))
  cat(sprintf("Outcome variable: %s\n", object$outcome_var))
  wt_str <- if (is.null(object$weight_var)) {
    "None (equal weights)"
  } else {
    object$weight_var
  }
  cat(sprintf("Weight variable: %s\n", wt_str))
  cat(sprintf("Frequency: %s\n", object$frequency))
  if (object$n_excluded_cells > 0L) {
    cat(sprintf("\nExcluded cells: %d\n", object$n_excluded_cells))
  }
  invisible(object)
}


#' @export
as.data.frame.lwdid_aggregation_result <- function(
  x, row.names = NULL, optional = FALSE, ...
) {
  as.data.frame(x$panel_data, row.names = row.names,
                optional = optional)
}


#' @export
to_csv.lwdid_aggregation_result <- function(
  x, path, include_metadata = TRUE, ...
) {
  if (include_metadata) {
    lines <- c(
      "# Aggregation Metadata",
      sprintf("# Original observations: %d", x$n_original_obs),
      sprintf("# Output cells: %d", x$n_cells),
      sprintf("# Units: %d", x$n_units),
      sprintf("# Periods: %d", x$n_periods),
      sprintf("# Unit variable: %s", x$unit_var),
      sprintf("# Time variable: %s",
              paste(x$time_var, collapse = ", ")),
      sprintf("# Outcome variable: %s", x$outcome_var),
      sprintf("# Weight variable: %s",
              if (is.null(x$weight_var)) "None" else x$weight_var),
      sprintf("# Frequency: %s", x$frequency),
      "#"
    )
    writeLines(lines, path)
    data.table::fwrite(x$panel_data, path, append = TRUE,
                        col.names = TRUE)
  } else {
    data.table::fwrite(x$panel_data, path)
  }
  invisible(x)
}


#' @export
to_dict.lwdid_aggregation_result <- function(x, ...) {
  list(
    n_original_obs = x$n_original_obs,
    n_cells = x$n_cells,
    n_units = x$n_units,
    n_periods = x$n_periods,
    min_cell_size = x$min_cell_size_stat,
    max_cell_size = x$max_cell_size,
    mean_cell_size = x$mean_cell_size,
    median_cell_size = x$median_cell_size,
    unit_var = x$unit_var,
    time_var = x$time_var,
    outcome_var = x$outcome_var,
    weight_var = x$weight_var,
    frequency = x$frequency,
    n_excluded_cells = x$n_excluded_cells
  )
}


# ============================================================================
# Story 11.3: Main Aggregation Function & Helpers
# ============================================================================

#' Build groupby columns based on frequency
#' @keywords internal
.build_group_columns <- function(unit_var, time_var, frequency) {
  c(unit_var, time_var)
}


#' Count unique periods in panel data
#' @keywords internal
.count_periods <- function(panel_data, time_var, frequency) {
  if (length(time_var) == 1L) {
    data.table::uniqueN(panel_data[[time_var]])
  } else {
    nrow(unique(panel_data[, ..time_var]))
  }
}


#' Build a single row for the aggregated panel
#' @param group_key Named list with group column values.
#' @param group_cols Character vector of grouping column names.
#' @param mean_val Numeric, weighted mean of the outcome.
#' @param n_obs Integer, number of observations in the cell.
#' @param group_data data.frame, raw data for this cell.
#' @param outcome_var Character, outcome variable name.
#' @param weight_var Character or NULL, weight variable name.
#' @param controls Character vector or NULL, control variable names.
#' @param treatment_var Character or NULL, treatment variable name.
#' @param gvar Character or NULL, cohort variable name.
#' @param variance Numeric or NULL, within-cell variance.
#' @param ess Numeric or NULL, effective sample size.
#' @param compute_variance Logical, whether to include variance column.
#' @param has_weights Logical, whether survey weights are used.
#' @return A named list representing one row of the aggregated panel.
#' @keywords internal
.build_aggregated_row <- function(group_key, group_cols, mean_val,
                                   n_obs, group_data, outcome_var,
                                   weight_var, controls,
                                   treatment_var, gvar,
                                   variance, ess,
                                   compute_variance, has_weights) {
  row <- list()

  # Group column values
  for (col in group_cols) {
    row[[col]] <- group_key[[col]]
  }

  # Outcome mean
  row[[outcome_var]] <- mean_val

  # Cell size
  row[["_n_obs"]] <- n_obs

  # Variance (always include column when compute_variance=TRUE)
  if (compute_variance) {
    row[[paste0(outcome_var, "_var")]] <- if (!is.null(variance)) variance else NA_real_
  }

  # ESS (always include column when survey weights used)
  if (has_weights) {
    row[["_ess"]] <- if (!is.null(ess)) ess else NA_real_
  }

  # Aggregate control variables
  if (!is.null(controls)) {
    raw_w <- if (!is.null(weight_var)) {
      group_data[[weight_var]]
    } else {
      NULL
    }
    for (ctrl in controls) {
      row[[ctrl]] <- .compute_control_weighted_mean(
        group_data[[ctrl]], raw_w
      )
    }
  }

  # Treatment variable (constant within cell)
  if (!is.null(treatment_var)) {
    row[[treatment_var]] <- group_data[[treatment_var]][1L]
  }

  # gvar (constant within unit)
  if (!is.null(gvar)) {
    row[[gvar]] <- group_data[[gvar]][1L]
  }

  row
}


#' Aggregate repeated cross-sectional data to panel format
#'
#' @description
#' Aggregates lower-level repeated cross-sectional data to the
#' unit-by-period level using weighted means, following
#' Lee & Wooldridge (2026) Section 3.
#'
#' @param data data.frame or data.table, repeated cross-sectional data
#' @param unit_var character, aggregation unit column name (e.g., "state")
#' @param time_var character or character vector, time variable column name(s)
#' @param outcome_var character, outcome variable column name
#' @param weight_var character or NULL, survey weight column name
#' @param controls character vector or NULL, control variable column names
#' @param treatment_var character or NULL, treatment indicator column name
#' @param gvar character or NULL, treatment timing variable column name
#' @param frequency character, aggregation frequency
#' @param min_cell_size positive integer, minimum observations per cell
#' @param compute_variance logical, whether to compute within-cell variance
#'
#' @return lwdid_aggregation_result S3 object
#' @export
aggregate_to_panel <- function(
  data,
  unit_var,
  time_var,
  outcome_var,
  weight_var = NULL,
  controls = NULL,
  treatment_var = NULL,
  gvar = NULL,
  frequency = c("annual", "quarterly", "monthly", "weekly"),
  min_cell_size = 1L,
  compute_variance = FALSE
) {
  # Step 1: Input validation
  .validate_aggregation_inputs(
    data, unit_var, time_var, outcome_var,
    weight_var, controls
  )

  # Validate frequency
  frequency_input <- tolower(frequency)
  frequency <- tryCatch(
    match.arg(frequency_input,
              c("annual", "quarterly", "monthly", "weekly")),
    error = function(e) {
      stop_lwdid(
        sprintf("Invalid frequency '%s'. Must be one of: annual, quarterly, monthly, weekly",
                frequency_input),
        class = "lwdid_invalid_parameter",
        param = "frequency",
        value = frequency_input,
        allowed = c("annual", "quarterly", "monthly", "weekly")
      )
    }
  )

  # Validate min_cell_size
  if (!is.numeric(min_cell_size) ||
      length(min_cell_size) != 1L ||
      min_cell_size < 1 ||
      min_cell_size != round(min_cell_size)) {
    stop_lwdid(
      sprintf("min_cell_size must be a positive integer. Got: %s",
              deparse(min_cell_size)),
      class = "lwdid_invalid_parameter",
      param = "min_cell_size",
      value = min_cell_size
    )
  }
  min_cell_size <- as.integer(min_cell_size)

  # Step 2: Data preparation
  data_work <- data.table::as.data.table(data.table::copy(data))
  n_original <- nrow(data_work)

  if (!is.null(weight_var)) {
    wt_result <- .validate_weights(data_work, weight_var)
    data_work <- wt_result$data
  }

  # Step 3: Treatment consistency validation
  if (!is.null(treatment_var)) {
    .validate_treatment_consistency(
      data_work, unit_var, time_var, treatment_var
    )
  }

  # Step 4: gvar consistency validation
  if (!is.null(gvar)) {
    .validate_gvar_consistency(data_work, unit_var, gvar)
  }

  # Step 5: Build group columns
  group_cols <- .build_group_columns(unit_var, time_var, frequency)

  # Step 6: Grouped aggregation
  unique_groups <- unique(data_work[, ..group_cols])
  n_groups <- nrow(unique_groups)
  aggregated_rows <- vector("list", n_groups)
  cell_stats_list <- vector("list", n_groups)
  excluded_cells <- list()
  row_idx <- 0L

  for (i in seq_len(n_groups)) {
    group_key <- as.list(unique_groups[i])
    group_data <- data_work[unique_groups[i],
                            on = group_cols, nomatch = NULL]

    # Compute weighted mean
    outcomes <- group_data[[outcome_var]]
    raw_w <- if (!is.null(weight_var)) {
      group_data[[weight_var]]
    } else {
      NULL
    }
    cell_result <- .compute_cell_weighted_mean(
      outcomes, raw_w, compute_variance
    )

    mean_val <- cell_result$mean
    variance <- cell_result$variance
    ess <- cell_result$ess
    n_obs <- cell_result$n_obs

    # Check cell size
    if (n_obs < min_cell_size) {
      excluded_cells <- c(excluded_cells, list(list(
        cell = group_key,
        n_obs = n_obs,
        reason = sprintf("below min_cell_size (%d)", min_cell_size)
      )))
      next
    }

    # Check all-NA/NaN outcome
    if (is.na(mean_val)) {
      excluded_cells <- c(excluded_cells, list(list(
        cell = group_key,
        n_obs = n_obs,
        reason = "all-NaN outcome"
      )))
      next
    }

    # Build output row
    row <- .build_aggregated_row(
      group_key, group_cols, mean_val, n_obs,
      group_data, outcome_var, weight_var, controls,
      treatment_var, gvar, variance, ess,
      compute_variance = compute_variance,
      has_weights = !is.null(weight_var)
    )
    row_idx <- row_idx + 1L
    aggregated_rows[[row_idx]] <- row

    # Build CellStatistics
    cell_stats_list[[row_idx]] <- new_cell_statistics(
      unit = group_key[[unit_var]],
      period = if (length(time_var) == 1L) {
        group_key[[time_var]]
      } else {
        group_key[time_var]
      },
      n_obs = n_obs,
      outcome_mean = mean_val,
      outcome_variance = variance,
      effective_sample_size = ess,
      weight_type = if (!is.null(weight_var)) "survey" else "equal"
    )
  }

  # Truncate pre-allocated lists
  if (row_idx > 0L) {
    aggregated_rows <- aggregated_rows[seq_len(row_idx)]
    cell_stats_list <- cell_stats_list[seq_len(row_idx)]
  } else {
    aggregated_rows <- list()
    cell_stats_list <- list()
  }

  # Step 7: All-excluded check
  if (length(aggregated_rows) == 0L) {
    max_obs <- if (length(excluded_cells) > 0L) {
      max(vapply(excluded_cells, function(x) x$n_obs, integer(1)))
    } else {
      0L
    }
    stop_lwdid(
      sprintf(
        paste0("All cells have fewer than %d observations ",
               "or all-NaN outcomes. ",
               "Total cells attempted: %d. ",
               "Consider reducing min_cell_size parameter."),
        min_cell_size, length(excluded_cells)
      ),
      class = c("lwdid_insufficient_cell_size",
                 "lwdid_aggregation_error"),
      min_cell_size = min_cell_size,
      max_observed_size = max_obs,
      n_cells = length(excluded_cells)
    )
  }

  # Step 8: Excluded cells warning
  n_excluded <- length(excluded_cells)
  if (n_excluded > 0L) {
    n_below_min <- sum(vapply(
      excluded_cells,
      function(x) grepl("min_cell_size", x$reason),
      logical(1)
    ))
    n_all_nan <- sum(vapply(
      excluded_cells,
      function(x) grepl("NaN", x$reason),
      logical(1)
    ))
    warn_lwdid(
      sprintf(
        "Excluded %d cells: %d below min_cell_size, %d with all-NaN outcomes.",
        n_excluded, n_below_min, n_all_nan
      ),
      class = "lwdid_data",
      detail = "excluded_cells",
      action_taken = "cells_removed"
    )
  }

  # Step 9: Build output panel
  panel_data <- data.table::rbindlist(aggregated_rows)

  # Step 10: Build CellStatistics table
  cell_stats_df <- data.table::rbindlist(cell_stats_list)

  # Step 11: Summary statistics
  cell_sizes <- vapply(
    cell_stats_list,
    function(x) x$n_obs, integer(1)
  )
  n_units <- data.table::uniqueN(panel_data[[unit_var]])
  n_periods <- .count_periods(panel_data, time_var, frequency)

  # Step 12: Small sample warning
  if (n_units < 3L) {
    warn_lwdid(
      sprintf(
        paste0("Aggregation resulted in %d units. ",
               "lwdid requires at least 3 units ",
               "for valid estimation."),
        n_units
      ),
      class = "lwdid_small_sample",
      n = n_units
    )
  }

  # Step 13: Build and return result
  new_lwdid_aggregation_result(
    panel_data = panel_data,
    n_original_obs = n_original,
    n_cells = length(aggregated_rows),
    n_units = n_units,
    n_periods = n_periods,
    cell_stats = cell_stats_df,
    min_cell_size_stat = min(cell_sizes),
    max_cell_size = max(cell_sizes),
    mean_cell_size = mean(cell_sizes),
    median_cell_size = stats::median(cell_sizes),
    unit_var = unit_var,
    time_var = time_var,
    outcome_var = outcome_var,
    weight_var = weight_var,
    frequency = frequency,
    n_excluded_cells = n_excluded,
    excluded_cells_info = excluded_cells
  )
}
