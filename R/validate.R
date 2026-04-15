# ============================================================================
# validate.R — Input validation framework for lwdid
#
# This file implements the complete input validation pipeline for the lwdid()
# function. All validation logic is centralized in validate_inputs(), which
# internally dispatches to layer-specific helper functions.
#
# Validation layers (executed in strict order):
#   Layer 1: Reserved column names + empty data check
#   Layer 2: Parameter type/range validation
#   Layer 3a: Mode identification (common_timing vs staggered)
#   Layer 3b: Treatment indicator binarization (on data.table copy)
#   Layer 4: Outcome/controls data type validation
#   Layer 5: Time-invariance validation (CT mode only)
#   Layer 6: Cross-parameter consistency checks
#   Layer 7: Data format conversion (data.table + string ID encoding)
#   Layer 8: Missing value handling (dropna)
#   Layer 9: Time index creation + staggered data validation
#   Layer 10: Panel structure validation
#
# Dependencies: conditions.R (Story E1-02), warning_registry.R (Story E1-03)
# ============================================================================

# ============================================================================
# Constants for validation thresholds
# ============================================================================

TIME_INVARIANCE_THRESHOLD <- 1e-10

COHORT_FLOAT_TOLERANCE <- 1e-9

NEAR_ZERO_TOLERANCE <- 1e-10

RESERVED_COLUMNS <- c("d_", "post_", "tindex", "tq",
                       "ydot", "ydot_postavg", "firstpost")

VALID_ROLLING_METHODS <- c("demean", "detrend", "demeanq", "detrendq")

VALID_ESTIMATORS <- c("ra", "ipw", "ipwra", "psm")

VALID_VCE_TYPES <- c("hc0", "hc1", "hc2", "hc3", "hc4",
                      "robust", "cluster", "bootstrap")

VALID_CONTROL_GROUPS <- c("not_yet_treated", "never_treated",
                           "all_others", "auto")

VALID_AGGREGATES <- c("none", "cohort", "overall", "event_time")

VALID_BALANCED_PANEL <- c("warn", "error", "ignore")

VALID_VERBOSE <- c("quiet", "default", "verbose")

VALID_MATCH_ORDER <- c("data", "random", "largest", "smallest")

VALID_TRIM_METHODS <- c("clip", "drop")

VALID_CALIPER_SCALES <- c("sd", "absolute")

VALID_RI_METHODS <- c("bootstrap", "permutation")

FREQ_LABELS <- list("4" = "quarter", "12" = "month", "52" = "week")


# ============================================================================
# Layer 1: Reserved column names + empty data check
# ============================================================================

#' Check that input data does not contain reserved column names
#'
#' @param data data.frame or data.table to check
#' @return NULL (invisible). Raises lwdid_invalid_parameter on conflict.
#' @keywords internal
.validate_reserved_columns <- function(data) {
  existing <- intersect(RESERVED_COLUMNS, names(data))
  if (length(existing) > 0L) {
    stop_lwdid(
      message = paste0(
        "Input data contains reserved column names: ",
        paste(existing, collapse = ", "), ".\n",
        "These columns are used internally by lwdid and will be overwritten.\n",
        "Please rename these columns before calling lwdid().\n",
        "Reserved column names: ",
        paste(RESERVED_COLUMNS, collapse = ", ")
      ),
      class = "lwdid_invalid_parameter",
      param = "data",
      value = existing,
      allowed = "columns not in reserved list"
    )
  }
  invisible(NULL)
}

#' Check that input data is a data.frame/data.table and non-empty
#'
#' @param data Input data to check
#' @return NULL (invisible). Raises lwdid_invalid_parameter if invalid.
#' @keywords internal
.validate_data_not_empty <- function(data) {
  if (!is.data.frame(data)) {
    stop_lwdid(
      message = sprintf(
        "Input data must be a data.frame or data.table. Got: %s",
        class(data)[1]
      ),
      class = "lwdid_invalid_parameter",
      param = "data",
      value = class(data)[1],
      allowed = "data.frame or data.table"
    )
  }
  if (nrow(data) < 1L) {
    stop_lwdid(
      message = "Input data is empty (0 rows).",
      class = "lwdid_invalid_parameter",
      param = "data",
      value = "0 rows",
      allowed = "nrow(data) >= 1"
    )
  }
  invisible(NULL)
}

# ============================================================================
# Layer 2: Parameter type/range validation
# ============================================================================

#' Validate a required string parameter that must be a column name in data
#' @param value The parameter value to validate
#' @param name Character string, parameter name for error messages
#' @param data data.frame to check column existence
#' @keywords internal
.validate_string_param <- function(value, name, data) {
  if (!is.character(value) || length(value) != 1L) {
    stop_lwdid(
      message = sprintf("'%s' must be a single character string. Got: %s",
                        name, class(value)[1]),
      class = "lwdid_invalid_parameter",
      param = name, value = value, allowed = "character(1)"
    )
  }
  if (!value %in% names(data)) {
    stop_lwdid(
      message = sprintf(
        "Required column '%s' not found in data. Available columns: %s",
        value, paste(names(data), collapse = ", ")),
      class = "lwdid_missing_column",
      column = value, available = names(data)
    )
  }
}

#' Validate optional string parameter (NULL allowed, column must exist if non-NULL)
#' @keywords internal
.validate_optional_string_param <- function(value, name, data) {
  if (is.null(value)) return(invisible(NULL))
  .validate_string_param(value, name, data)
}

#' Validate tvar parameter (length 1 or 2 character vector)
#' @keywords internal
.validate_tvar <- function(tvar, data) {
  if (!is.character(tvar) || !(length(tvar) %in% c(1L, 2L))) {
    stop_lwdid(
      message = paste0(
        "'tvar' must be a character string (annual) or length-2 character ",
        "vector (quarterly c(year_var, quarter_var)). Got: length ",
        length(tvar), ", class ", class(tvar)[1]
      ),
      class = "lwdid_invalid_parameter",
      param = "tvar", value = tvar,
      allowed = "character(1) or character(2)"
    )
  }
  for (tv in tvar) {
    if (!tv %in% names(data)) {
      stop_lwdid(
        message = sprintf(
          "Time variable column '%s' not found in data. Available: %s",
          tv, paste(names(data), collapse = ", ")),
        class = "lwdid_missing_column",
        column = tv, available = names(data)
      )
    }
    if (!is.numeric(data[[tv]])) {
      stop_lwdid(
        message = sprintf("Time variable '%s' must be numeric. Got: %s",
                          tv, class(data[[tv]])[1]),
        class = "lwdid_invalid_parameter",
        param = "tvar", value = class(data[[tv]])[1],
        allowed = "numeric"
      )
    }
  }
  # Quarterly data: validate quarter values in {1,2,3,4}
  if (length(tvar) == 2L) {
    quarter_vals <- unique(data[[tvar[2]]])
    quarter_vals <- quarter_vals[!is.na(quarter_vals)]
    if (!all(quarter_vals %in% c(1, 2, 3, 4))) {
      stop_lwdid(
        message = sprintf(
          "Quarter variable '%s' contains invalid values. Must be in {1,2,3,4}. Found: %s",
          tvar[2], paste(sort(quarter_vals), collapse = ", ")),
        class = "lwdid_invalid_parameter",
        param = "tvar", value = sort(quarter_vals),
        allowed = c(1, 2, 3, 4)
      )
    }
  }
}

#' Validate that a parameter value is in a set of valid choices
#' @keywords internal
.validate_choice <- function(value, valid_set, name) {
  if (!is.character(value) || length(value) != 1L ||
      !value %in% valid_set) {
    stop_lwdid(
      message = sprintf(
        "Invalid parameter '%s': got '%s'. Allowed values: %s",
        name, as.character(value),
        paste(valid_set, collapse = ", ")),
      class = "lwdid_invalid_parameter",
      param = name, value = value, allowed = valid_set
    )
  }
}

#' Validate VCE parameter (NULL or valid string)
#' @keywords internal
.validate_vce <- function(vce) {
  if (is.null(vce)) return(invisible(NULL))
  if (!is.character(vce) || length(vce) != 1L ||
      !vce %in% VALID_VCE_TYPES) {
    stop_lwdid(
      message = sprintf(
        "Invalid VCE type '%s'. Must be one of: NULL, %s",
        as.character(vce),
        paste(VALID_VCE_TYPES, collapse = ", ")),
      class = c("lwdid_invalid_vce", "lwdid_invalid_parameter"),
      vce_type = vce, allowed = c("NULL", VALID_VCE_TYPES)
    )
  }
}

#' Validate numeric in open/closed interval
#' @param exclusive logical. If TRUE, boundaries are exclusive (open interval).
#' @keywords internal
.validate_numeric_range <- function(value, lower, upper, name,
                                    exclusive = TRUE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    stop_lwdid(
      message = sprintf("'%s' must be a single numeric value. Got: %s",
                        name, class(value)[1]),
      class = "lwdid_invalid_parameter",
      param = name, value = value,
      allowed = sprintf("numeric in (%s, %s)", lower, upper)
    )
  }
  if (exclusive) {
    if (value <= lower || value >= upper) {
      stop_lwdid(
        message = sprintf(
          "'%s' must be in open interval (%s, %s). Got: %s",
          name, lower, upper, value),
        class = "lwdid_invalid_parameter",
        param = name, value = value,
        allowed = sprintf("(%s, %s)", lower, upper)
      )
    }
  } else {
    if (value < lower || value > upper) {
      stop_lwdid(
        message = sprintf(
          "'%s' must be in closed interval [%s, %s]. Got: %s",
          name, lower, upper, value),
        class = "lwdid_invalid_parameter",
        param = name, value = value,
        allowed = sprintf("[%s, %s]", lower, upper)
      )
    }
  }
}

#' Validate positive integer (>= min_val, default 1)
#' @keywords internal
.validate_positive_integer <- function(value, name, min_val = 1L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value != as.integer(value) || value < min_val) {
    stop_lwdid(
      message = sprintf(
        "'%s' must be an integer >= %d. Got: %s",
        name, min_val, as.character(value)),
      class = "lwdid_invalid_parameter",
      param = name, value = value,
      allowed = sprintf("integer >= %d", min_val)
    )
  }
}

#' Validate non-negative integer (>= 0)
#' @keywords internal
.validate_nonneg_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value != as.integer(value) || value < 0L) {
    stop_lwdid(
      message = sprintf(
        "'%s' must be a non-negative integer. Got: %s",
        name, as.character(value)),
      class = "lwdid_invalid_parameter",
      param = name, value = value,
      allowed = "integer >= 0"
    )
  }
}

#' Validate logical scalar
#' @keywords internal
.validate_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop_lwdid(
      message = sprintf("'%s' must be TRUE or FALSE. Got: %s",
                        name, as.character(value)),
      class = "lwdid_invalid_parameter",
      param = name, value = value, allowed = "TRUE or FALSE"
    )
  }
}

#' Validate optional positive numeric (NULL allowed)
#' @keywords internal
.validate_optional_positive_numeric <- function(value, name) {
  if (is.null(value)) return(invisible(NULL))
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value <= 0) {
    stop_lwdid(
      message = sprintf("'%s' must be a positive number or NULL. Got: %s",
                        name, as.character(value)),
      class = "lwdid_invalid_parameter",
      param = name, value = value, allowed = "positive numeric or NULL"
    )
  }
}

#' Validate optional integer (NULL allowed)
#' @keywords internal
.validate_optional_integer <- function(value, name) {
  if (is.null(value)) return(invisible(NULL))
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value != as.integer(value)) {
    stop_lwdid(
      message = sprintf("'%s' must be an integer or NULL. Got: %s",
                        name, as.character(value)),
      class = "lwdid_invalid_parameter",
      param = name, value = value, allowed = "integer or NULL"
    )
  }
}

#' Validate controls/ps_controls column existence
#' @keywords internal
.validate_controls_columns <- function(controls, data) {
  if (is.null(controls)) return(invisible(NULL))
  if (!is.character(controls)) {
    stop_lwdid(
      message = sprintf(
        "controls must be a character vector or NULL. Got: %s",
        class(controls)[1]),
      class = "lwdid_invalid_parameter",
      param = "controls", value = class(controls)[1],
      allowed = "character vector or NULL"
    )
  }
  missing_cols <- setdiff(controls, names(data))
  if (length(missing_cols) > 0L) {
    stop_lwdid(
      message = sprintf(
        "Control variable columns not found in data: %s. Available: %s",
        paste(missing_cols, collapse = ", "),
        paste(names(data), collapse = ", ")),
      class = "lwdid_missing_column",
      column = missing_cols, available = names(data)
    )
  }
}

#' Validate graph_options (NULL or named list)
#' @keywords internal
.validate_graph_options <- function(graph_options) {
  if (is.null(graph_options)) return(invisible(NULL))
  if (!is.list(graph_options)) {
    stop_lwdid(
      message = sprintf(
        "graph_options must be a named list or NULL. Got: %s",
        class(graph_options)[1]),
      class = "lwdid_invalid_parameter",
      param = "graph_options", value = class(graph_options)[1],
      allowed = "named list or NULL"
    )
  }
}

#' Validate season_var column values are integers in 1..Q range
#' @keywords internal
.validate_season_var_values <- function(data, season_var, Q) {
  if (is.null(season_var)) return(invisible(NULL))
  season_vals <- data[[season_var]]
  season_vals <- season_vals[!is.na(season_vals)]
  valid_range <- seq_len(Q)
  invalid_vals <- setdiff(unique(season_vals), valid_range)
  if (length(invalid_vals) > 0L) {
    stop_lwdid(
      message = sprintf(
        "season_var '%s' contains invalid values: %s. Must be integers in 1:%d.",
        season_var, paste(sort(invalid_vals), collapse = ", "), Q),
      class = "lwdid_invalid_parameter",
      param = "season_var",
      value = sort(invalid_vals),
      allowed = paste0("integers 1:", Q)
    )
  }
}

# ============================================================================
# Layer 2 (continued): Rolling parameter validation
# ============================================================================

#' Validate rolling parameter (case-insensitive) and seasonal requirements
#' @return character(1). Standardized rolling parameter in lowercase.
#' @keywords internal
.validate_rolling_parameter <- function(rolling, tvar, season_var) {
  if (!is.character(rolling) || length(rolling) != 1L) {
    stop_lwdid(
      message = sprintf(
        "'rolling' must be a single character string. Got: %s",
        class(rolling)[1]),
      class = c("lwdid_invalid_rolling", "lwdid_invalid_parameter"),
      method = rolling, allowed = VALID_ROLLING_METHODS
    )
  }
  rolling_lower <- tolower(rolling)
  if (!rolling_lower %in% VALID_ROLLING_METHODS) {
    stop_lwdid(
      message = sprintf(
        "Invalid rolling method '%s'. Must be one of: %s",
        rolling, paste(VALID_ROLLING_METHODS, collapse = ", ")),
      class = c("lwdid_invalid_rolling", "lwdid_invalid_parameter"),
      method = rolling, allowed = VALID_ROLLING_METHODS
    )
  }
  # Seasonal methods require either tvar=[year, quarter] or season_var
  if (rolling_lower %in% c("demeanq", "detrendq")) {
    has_tvar_list <- is.character(tvar) && length(tvar) == 2L
    has_season_var <- !is.null(season_var)
    if (!has_tvar_list && !has_season_var) {
      stop_lwdid(
        message = paste0(
          "rolling('", rolling_lower, "') requires either:\n",
          "  1. tvar = c(year_col, quarter_col) (length-2 vector), or\n",
          "  2. season_var parameter specifying the seasonal column.\n\n",
          "Example: lwdid(..., rolling='", rolling_lower,
          "', season_var='quarter', Q=4)"
        ),
        class = c("lwdid_invalid_rolling", "lwdid_invalid_parameter"),
        method = rolling_lower, allowed = VALID_ROLLING_METHODS
      )
    }
  }
  rolling_lower
}

# ============================================================================
# Layer 3a: Mode identification
# ============================================================================

#' Identify estimation mode: Common Timing or Staggered
#' @return character(1). "common_timing" or "staggered".
#' @keywords internal
.identify_mode <- function(d, post, gvar) {
  if (!is.null(gvar)) {
    if (!is.null(d) || !is.null(post)) {
      warn_lwdid(
        message = paste0(
          "Staggered mode (gvar specified): ",
          "d and post parameters are ignored."),
        class = "lwdid_data",
        detail = "gvar takes priority over d/post",
        action_taken = "d and post ignored"
      )
    }
    return("staggered")
  }
  if (!is.null(d) && !is.null(post)) {
    return("common_timing")
  }
  stop_lwdid(
    message = paste0(
      "Must specify either gvar (Staggered mode) or both d and post ",
      "(Common Timing mode). Got: gvar=NULL, d=",
      if (is.null(d)) "NULL" else d,
      ", post=",
      if (is.null(post)) "NULL" else post
    ),
    class = "lwdid_invalid_parameter",
    param = "mode", value = "incomplete",
    allowed = "gvar OR (d + post)"
  )
}

# ============================================================================
# Layer 3b: Binarization (executes after Layer 7, before Layer 8)
# ============================================================================

#' Binarize a numeric vector preserving NA
#'
#' Converts non-zero values to 1, zero to 0, NA stays NA.
#' Matches Python: (d_numeric != 0).where(d_numeric.notna()).astype('Int64')
#'
#' @param x numeric vector (may contain NA)
#' @return integer vector: 0L, 1L, or NA_integer_
#' @keywords internal
.binarize_with_na <- function(x) {
  x_numeric <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x_numeric), NA_integer_,
         ifelse(x_numeric != 0, 1L, 0L))
}


# ============================================================================
# Layer 4: Data type validation
# ============================================================================

#' Validate outcome variable is numeric (warn if logical)
#'
#' @param data data.frame
#' @param y character(1). Outcome variable name.
#' @return NULL (invisible).
#' @keywords internal
.validate_outcome_dtype <- function(data, y) {
  col <- data[[y]]
  if (is.logical(col)) {
    warn_lwdid(
      message = sprintf(
        "Outcome variable '%s' has logical type (TRUE/FALSE). This will be treated as numeric (1/0).", y
      ),
      class = "lwdid_data",
      detail = sprintf("'%s' is logical", y),
      action_taken = "treated as numeric 0/1"
    )
    return(invisible(NULL))
  }
  if (!is.numeric(col)) {
    stop_lwdid(
      message = sprintf(
        "Outcome variable '%s' must be numeric type. Found: %s",
        y, class(col)[1]
      ),
      class = "lwdid_invalid_parameter",
      param = y, value = class(col)[1], allowed = "numeric"
    )
  }
  invisible(NULL)
}

#' Validate all control variables are numeric
#'
#' @param data data.frame
#' @param controls character vector or NULL
#' @return NULL (invisible).
#' @keywords internal
.validate_controls_dtype <- function(data, controls) {
  if (is.null(controls)) return(invisible(NULL))
  non_numeric <- character(0)
  for (ctrl in controls) {
    if (!is.numeric(data[[ctrl]])) {
      non_numeric <- c(non_numeric,
                       sprintf("'%s' (%s)", ctrl, class(data[[ctrl]])[1]))
    }
  }
  if (length(non_numeric) > 0L) {
    stop_lwdid(
      message = sprintf(
        "Control variables must be numeric. Found non-numeric:\n  %s",
        paste(non_numeric, collapse = "\n  ")
      ),
      class = "lwdid_invalid_parameter",
      param = "controls", value = non_numeric, allowed = "all numeric"
    )
  }
  invisible(NULL)
}


# ============================================================================
# Layer 5: Time-invariance validation (Common Timing only)
# ============================================================================

#' Validate treatment indicator is time-invariant within each unit
#'
#' Uses within-unit standard deviation with threshold > 1e-10
#' (not strict equality, to handle floating-point precision).
#' Matches Python _validate_treatment_time_invariance().
#'
#' @param data data.frame
#' @param d character(1). Treatment indicator column name.
#' @param ivar character(1). Unit identifier column name.
#' @return NULL (invisible). Raises lwdid_invalid_parameter if time-varying.
#' @keywords internal
.validate_treatment_time_invariance <- function(data, d, ivar) {
  # Compute within-unit std for treatment indicator
  within_std <- tapply(data[[d]], data[[ivar]], sd, na.rm = TRUE)
  max_std <- max(within_std, na.rm = TRUE)

  if (!is.na(max_std) && max_std > TIME_INVARIANCE_THRESHOLD) {
    varying_units <- names(within_std)[
      which(within_std > TIME_INVARIANCE_THRESHOLD)
    ]
    n_varying <- length(varying_units)
    examples <- head(varying_units, 3)
    example_stds <- within_std[examples]
    example_details <- paste(
      sprintf("  Unit %s: std = %.6f", examples, example_stds),
      collapse = "\n"
    )
    stop_lwdid(
      message = paste0(
        "Treatment indicator '", d, "' must be time-invariant ",
        "(constant within each unit).\n",
        "Found ", n_varying,
        " units with time-varying treatment status:\n",
        example_details,
        if (n_varying > 3) "\n  ..." else "",
        "\n\nThis method requires unit-level D_i, not time-varying W_it.\n",
        "D_i should be constant across all periods for each unit."
      ),
      class = "lwdid_invalid_parameter",
      param = d, value = "time-varying",
      allowed = "time-invariant within unit"
    )
  }
  invisible(NULL)
}

#' Validate control variables are time-invariant within each unit
#'
#' Same threshold logic as treatment time-invariance (sd > 1e-10).
#' Matches Python _validate_time_invariant_controls().
#'
#' @param data data.frame
#' @param ivar character(1). Unit identifier column name.
#' @param controls character vector or NULL.
#' @return NULL (invisible). Raises lwdid_invalid_parameter if time-varying.
#' @keywords internal
.validate_time_invariant_controls <- function(data, ivar, controls) {
  if (is.null(controls)) return(invisible(NULL))
  time_varying <- list()
  for (ctrl in controls) {
    within_std <- tapply(data[[ctrl]], data[[ivar]], sd, na.rm = TRUE)
    max_std <- max(within_std, na.rm = TRUE)
    if (!is.na(max_std) && max_std > TIME_INVARIANCE_THRESHOLD) {
      n_varying <- sum(
        within_std > TIME_INVARIANCE_THRESHOLD, na.rm = TRUE
      )
      time_varying[[ctrl]] <- list(
        max_std = max_std, n_varying = n_varying
      )
    }
  }
  if (length(time_varying) > 0L) {
    details <- paste(
      sprintf(
        "  '%s': max within-unit std = %.6f, %d units time-varying",
        names(time_varying),
        vapply(time_varying, `[[`, numeric(1), "max_std"),
        vapply(time_varying, `[[`, numeric(1), "n_varying")
      ),
      collapse = "\n"
    )
    warn_lwdid(
      message = paste0(
        "Control variables are time-varying within units. ",
        "Pre-period means will be used as time-invariant proxy.\n",
        "Time-varying controls:\n", details
      ),
      class = "lwdid_data",
      detail = "controls_time_varying_detected",
      action_taken = "pre-period means will be used"
    )
  }
  invisible(NULL)
}


# ============================================================================
# Layer 6: Cross-parameter consistency checks (Table J rules)
# ============================================================================

#' Validate cross-parameter consistency (Table J rules)
#'
#' Checks logical dependencies between parameters.
#' Must be called after individual parameter validation (Layer 2).
#'
#' @param vce,cluster_var,rolling,season_var,tvar,estimator Parameters
#' @param controls,ps_controls,aggregate,control_group Parameters
#' @param graph,ri,rireps,ri_method,exclude_pre_periods,Q,mode Parameters
#' @return NULL (invisible). Raises errors/warnings as needed.
#' @keywords internal
.validate_cross_param_consistency <- function(
  vce, cluster_var, rolling, season_var, tvar, estimator,
  controls, ps_controls, aggregate, control_group,
  graph, ri, rireps, ri_method, exclude_pre_periods, Q, mode
) {
  # J1: vce requires cluster_var
  if ((!is.null(vce)) && vce %in% c("cluster", "bootstrap") &&
      is.null(cluster_var)) {
    stop_lwdid(
      message = sprintf(
        "vce='%s' requires cluster_var to be specified.", vce
      ),
      class = "lwdid_invalid_parameter",
      param = "cluster_var", value = "NULL",
      allowed = "non-NULL when vce is 'cluster' or 'bootstrap'"
    )
  }

  # J3: ps_controls only effective for PS estimators
  if (estimator == "ra" && !is.null(ps_controls)) {
    warn_lwdid(
      message = paste0(
        "ps_controls specified but estimator='ra'. ",
        "ps_controls are ignored for RA estimator."
      ),
      class = "lwdid_data",
      detail = "ps_controls not used with estimator='ra'",
      action_taken = "ps_controls ignored"
    )
  }

  # J4: PSM needs controls
  if (estimator == "psm" &&
      is.null(controls) && is.null(ps_controls)) {
    stop_lwdid(
      message = paste0(
        "estimator='psm' requires controls or ",
        "ps_controls to be specified."
      ),
      class = "lwdid_invalid_parameter",
      param = "controls", value = "NULL",
      allowed = "non-NULL for PSM estimator"
    )
  }

  # J5: IPWRA needs controls
  if (estimator == "ipwra" && is.null(controls)) {
    stop_lwdid(
      message = paste0(
        "estimator='ipwra' requires controls (outcome model). ",
        "ps_controls cannot substitute."
      ),
      class = "lwdid_invalid_parameter",
      param = "controls", value = "NULL",
      allowed = "non-NULL for IPWRA estimator"
    )
  }

  # J6: IPW needs controls or ps_controls
  if (estimator == "ipw" &&
      is.null(controls) && is.null(ps_controls)) {
    stop_lwdid(
      message = paste0(
        "estimator='ipw' requires controls or ps_controls ",
        "for propensity score model."
      ),
      class = "lwdid_invalid_parameter",
      param = "controls", value = "NULL",
      allowed = "non-NULL for IPW estimator"
    )
  }

  # J10: graph requires ggplot2
  if (isTRUE(graph) &&
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop_lwdid(
      message = paste0(
        "graph=TRUE requires the ggplot2 package. ",
        "Install with: install.packages('ggplot2')"
      ),
      class = "lwdid_visualization_error",
      plot_type = "lwdid", detail = "ggplot2 not installed"
    )
  }

  # J11: ri_method/rireps only effective when ri=TRUE (silently ignored)

  # J12: cluster_var without cluster vce
  if ((!is.null(vce)) &&
      (!vce %in% c("cluster", "bootstrap")) &&
      !is.null(cluster_var)) {
    warn_lwdid(
      message = sprintf(
        "cluster_var specified but vce='%s' does not use clustering. cluster_var ignored.",
        vce
      ),
      class = "lwdid_data",
      detail = "cluster_var not needed",
      action_taken = "cluster_var ignored"
    )
  }
  if (is.null(vce) && !is.null(cluster_var)) {
    warn_lwdid(
      message = "cluster_var specified but vce=NULL. cluster_var ignored.",
      class = "lwdid_data",
      detail = "cluster_var not needed",
      action_taken = "cluster_var ignored"
    )
  }
}


# ============================================================================
# Layer 7: Data format conversion
# ============================================================================

#' Convert input data to data.table (copy semantics)
#'
#' If data is already a data.table, creates a copy.
#' If data is a data.frame, converts to data.table (implicit copy).
#' Never modifies the original data object.
#'
#' @param data data.frame or data.table
#' @return data.table copy
#' @keywords internal
.convert_to_datatable <- function(data) {
  if (data.table::is.data.table(data)) {
    data.table::copy(data)
  } else {
    data.table::as.data.table(data)
  }
}

#' Convert string unit identifiers to numeric codes
#'
#' If ivar column is character/factor, converts to consecutive integers
#' starting from 1. Returns bidirectional mapping for recovery.
#' Matches Python _convert_string_id().
#'
#' NOTE: R uses sort() for deterministic alphabetical ordering, while
#' Python pd.factorize() preserves order of appearance.
#'
#' @param dt data.table (modified in place)
#' @param ivar character(1). Unit identifier column name.
#' @return list with $data (data.table) and $id_mapping (list or NULL)
#' @keywords internal
.convert_string_id <- function(dt, ivar) {
  id_mapping <- NULL
  col <- dt[[ivar]]
  if (is.character(col) || is.factor(col)) {
    col_char <- as.character(col)
    uniques <- unique(col_char[!is.na(col_char)])
    uniques <- sort(uniques)  # deterministic ordering
    original_to_numeric <- setNames(seq_along(uniques), uniques)
    numeric_to_original <- setNames(uniques, seq_along(uniques))
    # Map values
    data.table::set(dt, j = ivar,
                    value = original_to_numeric[col_char])
    id_mapping <- list(
      original_to_numeric = as.list(original_to_numeric),
      numeric_to_original = as.list(numeric_to_original)
    )
  }
  list(data = dt, id_mapping = id_mapping)
}


# ============================================================================
# Layer 8: Missing value handling
# ============================================================================

#' Drop rows with NA in required variables and warn
#'
#' Required variables: y, ivar, tvar columns, d_/post_ (Common Timing).
#' Must be called BEFORE .create_time_index() (Layer 9).
#' Matches Python dropna + DataWarning.
#'
#' @param dt data.table
#' @param y character(1). Outcome variable name.
#' @param ivar character(1). Unit identifier column name.
#' @param tvar character(1) or character(2). Time variable(s).
#' @param mode character(1). "common_timing" or "staggered".
#' @param d character(1) or NULL. Treatment indicator column name.
#' @param post character(1) or NULL. Post indicator column name.
#' @return data.table with NA rows removed
#' @keywords internal
.handle_missing_values <- function(dt, y, ivar, tvar, mode, d, post) {
  required_vars <- c(y, ivar, tvar)
  if (mode == "common_timing") {
    required_vars <- c(required_vars, "d_", "post_")
  }

  n_before <- nrow(dt)
  na_mask <- rowSums(is.na(dt[, ..required_vars])) > 0L
  if (any(na_mask)) {
    dt <- dt[!na_mask]
    n_dropped <- sum(na_mask)
    # Build display names
    var_display <- required_vars
    if ("d_" %in% var_display) {
      var_display[var_display == "d_"] <- sprintf(
        "%s (treatment indicator)", d
      )
    }
    if ("post_" %in% var_display) {
      var_display[var_display == "post_"] <- sprintf(
        "%s (post indicator)", post
      )
    }
    warn_lwdid(
      message = sprintf(
        "Dropped %d observations due to missing values in required variables: %s",
        n_dropped, paste(var_display, collapse = ", ")
      ),
      class = "lwdid_data",
      detail = sprintf(
        "%d rows with NA in required variables", n_dropped
      ),
      action_taken = sprintf("dropped %d rows", n_dropped)
    )
  }
  dt
}


# ============================================================================
# Layer 9: Time index creation + staggered data validation
# ============================================================================

#' Create sequential time index (tindex) starting from 1
#'
#' For annual data: tindex = year - min(year) + 1
#' For quarterly data: tq = (year - 1960) * 4 + quarter;
#'   tindex = tq - min(tq) + 1
#' Must be called AFTER .handle_missing_values() (Layer 8).
#' Matches Python _create_time_index().
#'
#' @param dt data.table
#' @param tvar character(1) or character(2)
#' @return list with $data and $is_quarterly
#' @keywords internal
.create_time_index <- function(dt, tvar) {
  if (length(tvar) == 1L) {
    # Annual data
    tvar_name <- tvar
    tvar_vals <- dt[[tvar_name]]
    dt[, tindex := tvar_vals - min(tvar_vals, na.rm = TRUE) + 1L]
    is_quarterly <- FALSE
  } else {
    # Quarterly data
    year_name <- tvar[1]
    quarter_name <- tvar[2]
    dt[, tq := (dt[[year_name]] - 1960L) * 4L + dt[[quarter_name]]]
    dt[, tindex := tq - min(tq, na.rm = TRUE) + 1L]
    is_quarterly <- TRUE
  }
  list(data = dt, is_quarterly = is_quarterly)
}


# ============================================================================
# Temporary stub for is_never_treated() (Story E1-07)
# Remove this stub when Story E1-07 is implemented.
# ============================================================================

#' Validate staggered DiD data and extract cohort structure
#'
#' Performs comprehensive validation for staggered adoption settings.
#' Matches Python validate_staggered_data().
#'
#' @param dt data.table with tindex already created
#' @param gvar character(1). Cohort variable column name.
#' @param ivar character(1). Unit identifier column name.
#' @param tvar character(1) or character(2). Time variable(s).
#' @param y character(1). Outcome variable column name.
#' @return list with cohort structure information
#' @keywords internal
.validate_staggered_data <- function(dt, gvar, ivar, tvar, y) {
  # 1. gvar must be numeric
  if (!is.numeric(dt[[gvar]])) {
    stop_lwdid(
      message = sprintf(
        "gvar column '%s' must be numeric, got %s.",
        gvar, class(dt[[gvar]])[1]
      ),
      class = "lwdid_invalid_staggered_data",
      gvar = gvar, invalid_values = class(dt[[gvar]])[1],
      detail = "non-numeric gvar"
    )
  }

  # 2. gvar must be time-invariant within unit
  gvar_col <- gvar  # column name string
  gvar_nunique <- dt[, .(n = data.table::uniqueN(.SD[[1L]])),
                     by = c(ivar), .SDcols = gvar_col]
  inconsistent <- gvar_nunique[n > 1L]
  if (nrow(inconsistent) > 0L) {
    bad_units <- head(inconsistent[[ivar]], 5)
    stop_lwdid(
      message = paste0(
        "gvar must be time-invariant within each unit. ",
        "Units with varying gvar: ",
        paste(bad_units, collapse = ", "),
        if (nrow(inconsistent) > 5) "..." else ""
      ),
      class = "lwdid_invalid_staggered_data",
      gvar = gvar, invalid_values = bad_units,
      detail = "gvar varies within unit"
    )
  }

  # 3. Get unit-level gvar
  unit_gvar <- dt[, .(gval = .SD[[1L]][1L]), by = c(ivar), .SDcols = gvar_col]

  # 4. Check for negative values (including -Inf)
  non_na_gvar <- unit_gvar$gval[!is.na(unit_gvar$gval)]
  negative_gvar <- non_na_gvar[non_na_gvar < 0]
  if (length(negative_gvar) > 0L) {
    neg_vals <- head(sort(unique(negative_gvar)), 5)
    stop_lwdid(
      message = paste0(
        "gvar column contains negative values: ",
        paste(neg_vals, collapse = ", "), ". ",
        "Valid: positive integer (cohort), 0 (NT), ",
        "Inf (NT), NA (NT)."
      ),
      class = "lwdid_invalid_staggered_data",
      gvar = gvar, invalid_values = neg_vals,
      detail = "negative gvar values"
    )
  }

  # 5. Identify never-treated units
  unit_gvar[, is_nt := vapply(gval, is_never_treated, logical(1))]
  n_nt <- sum(unit_gvar$is_nt)
  treated_gvar <- unit_gvar[is_nt == FALSE]

  # 6. Must have at least one treatment cohort
  if (nrow(treated_gvar) == 0L) {
    stop_lwdid(
      message = paste0(
        "No treatment cohorts found. ",
        "All units are never-treated."
      ),
      class = "lwdid_invalid_staggered_data",
      gvar = gvar,
      invalid_values = head(unit_gvar$gval, 10),
      detail = "no valid treatment cohorts"
    )
  }

  # 7. Extract cohorts and sizes
  cohorts <- sort(unique(as.integer(treated_gvar$gval)))
  cohort_sizes <- as.list(table(as.integer(treated_gvar$gval)))
  n_treated <- nrow(treated_gvar)

  # 8. Time range
  tvar_col <- if (length(tvar) == 1L) tvar else tvar[1]
  T_min <- as.integer(min(dt[[tvar_col]], na.rm = TRUE))
  T_max <- as.integer(max(dt[[tvar_col]], na.rm = TRUE))

  # 9. Build warning list
  warning_list <- character(0)

  if (n_nt == 0L) {
    warning_list <- c(warning_list,
      paste0(
        "No never-treated units found. ",
        "Only (g,r)-specific effects can be estimated ",
        "using not-yet-treated controls. ",
        "Use aggregate='none'."
      )
    )
  }

  cohorts_outside <- cohorts[cohorts < T_min | cohorts > T_max]
  if (length(cohorts_outside) > 0L) {
    warning_list <- c(warning_list,
      sprintf(
        "Cohorts outside observed time range [%d, %d]: %s.",
        T_min, T_max,
        paste(cohorts_outside, collapse = ", ")
      )
    )
  }

  min_cohort <- min(cohorts)
  if (min_cohort <= T_min) {
    warning_list <- c(warning_list,
      sprintf(
        "Earliest cohort (%d) has no pre-treatment period (T_min=%d).",
        min_cohort, T_min
      )
    )
  }

  if (n_nt == 1L) {
    warning_list <- c(warning_list,
      paste0(
        "Only 1 never-treated unit. ",
        "Inference for cohort/overall effects may be unreliable."
      )
    )
  }

  # 10. Outcome missing values
  n_missing_y <- sum(is.na(dt[[y]]))
  if (n_missing_y > 0L) {
    pct <- n_missing_y / nrow(dt) * 100
    warning_list <- c(warning_list,
      sprintf(
        "Outcome '%s' has %d missing values (%.1f%%).",
        y, n_missing_y, pct
      )
    )
  }

  list(
    cohorts = cohorts,
    n_cohorts = length(cohorts),
    n_never_treated = n_nt,
    n_treated = n_treated,
    has_never_treated = n_nt > 0L,
    cohort_sizes = cohort_sizes,
    T_min = T_min,
    T_max = T_max,
    N_total = nrow(unit_gvar),
    N_obs = nrow(dt),
    warning_list = warning_list
  )
}


# ============================================================================
# Layer 10: Panel structure validation
# ============================================================================

#' Check for duplicate (ivar, tvar) observations
#'
#' @param dt data.table
#' @param ivar character(1)
#' @param tvar character(1) or character(2)
#' @return NULL (invisible). Raises lwdid_invalid_parameter if duplicates.
#' @keywords internal
.validate_no_duplicates <- function(dt, ivar, tvar) {
  if (length(tvar) == 2L) {
    dup_cols <- c(ivar, tvar[1], tvar[2])
  } else {
    dup_cols <- c(ivar, "tindex")
  }
  dup_mask <- duplicated(dt, by = dup_cols) |
    duplicated(dt, by = dup_cols, fromLast = TRUE)
  if (any(dup_mask)) {
    n_dup <- sum(dup_mask)
    stop_lwdid(
      message = sprintf(
        "Duplicate observations found. Each (unit, time) must be unique. Found %d duplicate rows.",
        n_dup
      ),
      class = "lwdid_invalid_parameter",
      param = "data",
      value = sprintf("%d duplicates", n_dup),
      allowed = "unique (ivar, tvar) combinations"
    )
  }
}


#' Validate Common Timing specific panel structure
#'
#' Checks: post cross-unit consistency, post monotonicity,
#' pre/post period existence, treated/control existence.
#'
#' @param dt data.table with tindex, d_, post_ columns
#' @param ivar character(1)
#' @return list with K, tpost1, n_treated, n_control
#' @keywords internal
.validate_common_timing_structure <- function(dt, ivar) {
  # Pre/post existence
  pre_obs <- dt[post_ == 0L]
  post_obs <- dt[post_ == 1L]
  if (nrow(pre_obs) == 0L) {
    stop_lwdid(
      message = "No pre-treatment observations (post==0).",
      class = "lwdid_insufficient_data",
      n = 0L, n_treated = 0L, n_control = 0L
    )
  }
  if (nrow(post_obs) == 0L) {
    stop_lwdid(
      message = "No post-treatment observations (post==1).",
      class = "lwdid_insufficient_data",
      n = 0L, n_treated = 0L, n_control = 0L
    )
  }

  K <- as.integer(max(dt[post_ == 0L, tindex]))
  tpost1 <- as.integer(min(dt[post_ == 1L, tindex]))

  # Post cross-unit consistency (common timing assumption)
  post_by_time <- dt[, .(n_unique = data.table::uniqueN(post_)),
                     by = tindex]
  violating <- post_by_time[n_unique > 1L]
  if (nrow(violating) > 0L) {
    stop_lwdid(
      message = sprintf(
        "Common timing assumption violated: post varies across units at tindex=%s.",
        paste(violating$tindex, collapse = ", ")
      ),
      class = "lwdid_invalid_parameter",
      param = "post", value = "varies across units",
      allowed = "constant across units at each time"
    )
  }

  # Post monotonicity
  post_by_time_val <- dt[, .(post_val = post_[1]),
                         by = tindex][order(tindex)]
  first_post_idx <- which(post_by_time_val$post_val == 1L)
  if (length(first_post_idx) > 0L) {
    first_post_t <- post_by_time_val$tindex[first_post_idx[1]]
    after_first <- post_by_time_val[tindex > first_post_t]
    reversals <- after_first[post_val == 0L]
    if (nrow(reversals) > 0L) {
      stop_lwdid(
        message = sprintf(
          "Post variable is not monotone: found post=0 at tindex=%s after first post=1 at tindex=%d.",
          paste(reversals$tindex, collapse = ", "),
          first_post_t
        ),
        class = "lwdid_time_discontinuity",
        gaps = reversals$tindex, ivar = ivar,
        tvar = "tindex"
      )
    }
  }

  # Treated/control counts
  unit_d <- dt[, .(d_max = max(d_, na.rm = TRUE)), by = c(ivar)]
  n_treated <- sum(unit_d$d_max == 1L)
  n_control <- sum(unit_d$d_max == 0L)

  if (n_treated < 1L) {
    stop_lwdid(
      message = "No treated units found (d==1).",
      class = c("lwdid_no_treated", "lwdid_insufficient_data"),
      treat_var = "d_", n_units = nrow(unit_d)
    )
  }
  if (n_control < 1L) {
    stop_lwdid(
      message = "No control units found (d==0).",
      class = c("lwdid_no_control", "lwdid_insufficient_data"),
      treat_var = "d_", n_units = nrow(unit_d)
    )
  }

  # Single treated/control warnings
  if (n_treated == 1L) {
    warn_lwdid(
      message = sprintf(
        "Only 1 treated unit (N_treated=1). Results may be unreliable."
      ),
      class = "lwdid_small_sample",
      n = nrow(unit_d),
      n_treated = n_treated, n_control = n_control
    )
  }
  if (n_control == 1L) {
    warn_lwdid(
      message = sprintf(
        "Only 1 control unit (N_control=1). Results may be unreliable."
      ),
      class = "lwdid_small_sample",
      n = nrow(unit_d),
      n_treated = n_treated, n_control = n_control
    )
  }

  list(K = K, tpost1 = tpost1,
       n_treated = n_treated, n_control = n_control)
}


#' Check panel balance and handle according to balanced_panel parameter
#'
#' @param dt data.table
#' @param ivar character(1)
#' @param balanced_panel character(1). "warn", "error", or "ignore".
#' @param mode character(1). "common_timing" or "staggered".
#' @param gvar character(1) or NULL.
#' @param tvar character(1) or character(2). Time variable(s).
#' @return logical. TRUE if balanced, FALSE if unbalanced.
#' @keywords internal
.validate_panel_balance <- function(dt, ivar, balanced_panel,
                                    mode, gvar, tvar) {
  panel_counts <- dt[, .N, by = c(ivar)]
  if (data.table::uniqueN(panel_counts$N) == 1L) {
    return(TRUE)  # Balanced
  }

  min_obs <- min(panel_counts$N)
  max_obs <- max(panel_counts$N)
  n_incomplete <- sum(panel_counts$N < max_obs)
  pct_unbalanced <- n_incomplete / nrow(panel_counts) * 100

  if (balanced_panel == "error") {
    stop_lwdid(
      message = sprintf(
        "Unbalanced panel: obs per unit ranges from %d to %d. %d units (%.1f%%) incomplete.",
        min_obs, max_obs, n_incomplete, pct_unbalanced
      ),
      class = "lwdid_unbalanced_panel",
      min_obs = min_obs, max_obs = max_obs,
      n_incomplete_units = n_incomplete,
      pct_unbalanced = pct_unbalanced
    )
  }

  if (balanced_panel == "warn") {
    msg <- sprintf(
      "Unbalanced panel: obs per unit ranges from %d to %d. %d units (%.1f%%) incomplete.",
      min_obs, max_obs, n_incomplete, pct_unbalanced
    )
    # Staggered mode: add method usability diagnostics
    if (mode == "staggered" && !is.null(gvar)) {
      usability <- .compute_method_usability(
        dt, ivar, gvar, tvar
      )
      msg <- paste0(msg, sprintf(
        "\nDiagnostics: %.1f%% usable for demean, %.1f%% usable for detrend.",
        usability$pct_demean, usability$pct_detrend
      ))
    }
    warn_lwdid(
      message = msg,
      class = "lwdid_data",
      detail = sprintf(
        "unbalanced panel: %d-%d obs", min_obs, max_obs
      ),
      action_taken = "continued with unbalanced panel"
    )
  }

  # balanced_panel == "ignore": silent
  return(FALSE)
}


#' Compute percentage of treated units usable for demean/detrend
#'
#' Uses original time variable (tvar) rather than tindex for comparison
#' with g, because g is in the original time scale (e.g. 2005) while
#' tindex starts from 1.
#' Matches Python _compute_method_usability().
#'
#' @param dt data.table
#' @param ivar character(1)
#' @param gvar character(1)
#' @param tvar character(1) or character(2). Time variable(s).
#' @return list with pct_demean and pct_detrend
#' @keywords internal
.compute_method_usability <- function(dt, ivar, gvar, tvar) {
  unit_gvar <- dt[, .(g = .SD[[1L]][1L]), by = c(ivar), .SDcols = gvar]
  unit_gvar[, is_nt := vapply(g, is_never_treated, logical(1))]
  treated <- unit_gvar[is_nt == FALSE]

  if (nrow(treated) == 0L) {
    return(list(pct_demean = 100, pct_detrend = 100))
  }

  n_treated <- nrow(treated)
  n_below_demean <- 0L
  n_below_detrend <- 0L

  # Use original time variable (not tindex) vs g
  tvar_col <- if (length(tvar) == 1L) tvar else tvar[1]
  for (i in seq_len(nrow(treated))) {
    uid <- treated[[ivar]][i]
    g <- treated$g[i]
    unit_data <- dt[dt[[ivar]] == uid]
    n_pre <- sum(unit_data[[tvar_col]] < g, na.rm = TRUE)
    if (n_pre < 1L) n_below_demean <- n_below_demean + 1L
    if (n_pre < 2L) n_below_detrend <- n_below_detrend + 1L
  }

  list(
    pct_demean = 100 * (1 - n_below_demean / n_treated),
    pct_detrend = 100 * (1 - n_below_detrend / n_treated)
  )
}


#' Validate pre-treatment period count meets rolling method requirements
#'
#' Requirements (from lw2025/lw2026):
#'   demean:   >= 1 (Procedure 2.1/4.1)
#'   detrend:  >= 2 (Procedure 3.1/5.1)
#'   demeanq:  >= Q+1 (Q params, need df >= 1)
#'   detrendq: >= Q+2 (Q+1 params, need df >= 1)
#'
#' @param mode character(1)
#' @param rolling character(1)
#' @param n_pre_periods integer(1). Number of pre-treatment periods.
#' @param exclude_pre_periods integer(1). Periods to exclude.
#' @param Q integer(1). Number of seasons.
#' @param cohorts integer vector. Treatment cohorts (Staggered).
#' @param T_min integer(1). Minimum time value.
#' @param dt data.table
#' @param ivar character(1)
#' @param gvar character(1) or NULL
#' @return NULL (invisible).
#' @keywords internal
.validate_pre_period_sufficiency <- function(
  mode, rolling, n_pre_periods, exclude_pre_periods, Q,
  cohorts, T_min, dt, ivar, gvar
) {
  min_required <- switch(rolling,
    demean = 1L,
    detrend = 2L,
    demeanq = Q + 1L,
    detrendq = Q + 2L,
    1L  # fallback
  )

  if (mode == "common_timing") {
    effective_pre <- n_pre_periods - exclude_pre_periods
    if (effective_pre < min_required) {
      stop_lwdid(
        message = sprintf(
          "Insufficient pre-treatment periods for '%s': available=%d, required=%d.%s",
          rolling, effective_pre, min_required,
          if (exclude_pre_periods > 0L)
            sprintf(" (excluded %d periods)", exclude_pre_periods)
          else ""
        ),
        class = c("lwdid_insufficient_pre_periods",
                   "lwdid_insufficient_data"),
        available = effective_pre,
        required = min_required,
        rolling = rolling, cohort = NULL,
        excluded = exclude_pre_periods
      )
    }
  } else {
    # Staggered: check per-cohort
    insufficient_cohorts <- list()
    for (g in cohorts) {
      n_pre_g <- g - T_min
      effective_pre_g <- n_pre_g - exclude_pre_periods
      if (effective_pre_g < min_required) {
        insufficient_cohorts[[as.character(g)]] <- effective_pre_g
      }
    }

    if (length(insufficient_cohorts) == length(cohorts)) {
      # ALL cohorts insufficient
      stop_lwdid(
        message = paste0(
          "All cohorts have insufficient pre-treatment periods for '",
          rolling, "'. Required: ", min_required, ". Cohorts: ",
          paste(sprintf("g=%s (available=%d)",
                        names(insufficient_cohorts),
                        unlist(insufficient_cohorts)),
                collapse = ", ")
        ),
        class = c("lwdid_insufficient_pre_periods",
                   "lwdid_insufficient_data"),
        available = 0L, required = min_required,
        rolling = rolling
      )
    } else if (length(insufficient_cohorts) > 0L) {
      # Some cohorts insufficient - warning
      warn_lwdid(
        message = paste0(
          "Some cohorts have insufficient pre-periods for '",
          rolling, "' (required=", min_required, "): ",
          paste(sprintf("g=%s (%d)",
                        names(insufficient_cohorts),
                        unlist(insufficient_cohorts)),
                collapse = ", "),
          ". These cohorts will be skipped."
        ),
        class = "lwdid_data",
        detail = "insufficient pre-periods for some cohorts",
        action_taken = "cohorts skipped"
      )
    }
  }
}


#' Validate seasonal diversity for demeanq/detrendq
#'
#' Checks: (1) each unit has >= 2 distinct seasons in pre-period,
#' (2) all post-period seasons appear in pre-period.
#' Supports both Common Timing and Staggered modes.
#'
#' @param dt data.table
#' @param ivar character(1)
#' @param season_var character(1). Season column name.
#' @param post_var character(1) or NULL. "post_" for CT; NULL for Staggered.
#' @param Q integer(1). Number of seasons per cycle.
#' @param gvar character(1) or NULL. Cohort variable (Staggered only).
#' @param tvar character or NULL. Time variable(s) (Staggered only).
#' @return NULL (invisible).
#' @keywords internal
.validate_season_diversity <- function(dt, ivar, season_var, post_var,
                                       Q, gvar = NULL, tvar = NULL) {
  freq_label <- FREQ_LABELS[[as.character(Q)]]
  if (is.null(freq_label)) freq_label <- "season"

  if (!is.null(post_var)) {
    # --- Common Timing mode: use post_var column ---
    units <- unique(dt[[ivar]])
    for (uid in units) {
      unit_pre <- dt[dt[[ivar]] == uid & dt[[post_var]] == 0L]
      unit_post <- dt[dt[[ivar]] == uid & dt[[post_var]] == 1L]
      .check_unit_season_diversity(
        uid, unit_pre, unit_post, season_var, freq_label
      )
    }
  } else if (!is.null(gvar) && !is.null(tvar)) {
    # --- Staggered mode ---
    tvar_col <- if (length(tvar) == 1L) tvar else tvar[1]
    unit_gvar <- dt[, .(g = .SD[[1L]][1L]), by = c(ivar), .SDcols = gvar]
    # Only check treated units (NT units skip)
    treated <- unit_gvar[
      !vapply(g, is_never_treated, logical(1))
    ]
    for (i in seq_len(nrow(treated))) {
      uid <- treated[[ivar]][i]
      g <- treated$g[i]
      unit_data <- dt[dt[[ivar]] == uid]
      unit_pre <- unit_data[unit_data[[tvar_col]] < g]
      unit_post <- unit_data[unit_data[[tvar_col]] >= g]
      .check_unit_season_diversity(
        uid, unit_pre, unit_post, season_var, freq_label
      )
    }
  }
}


#' Check season diversity for a single unit (shared helper)
#' @param uid Unit identifier value
#' @param unit_pre data.table of pre-period rows
#' @param unit_post data.table of post-period rows
#' @param season_var character(1). Season column name.
#' @param freq_label character(1). Label for error messages.
#' @keywords internal
.check_unit_season_diversity <- function(uid, unit_pre, unit_post,
                                          season_var, freq_label) {
  pre_seasons <- unique(unit_pre[[season_var]])
  pre_seasons <- pre_seasons[!is.na(pre_seasons)]

  if (length(pre_seasons) < 2L) {
    stop_lwdid(
      message = sprintf(
        "Unit %s has only %d %s(s) in pre-period. demeanq/detrendq requires >= 2 different %ss. Found: %s",
        as.character(uid), length(pre_seasons), freq_label,
        freq_label, paste(sort(pre_seasons), collapse = ", ")
      ),
      class = c("lwdid_insufficient_quarter_diversity",
                 "lwdid_insufficient_data"),
      unit = uid, missing_quarters = integer(0),
      pre_quarters = sort(pre_seasons),
      post_quarters = integer(0)
    )
  }

  post_seasons <- unique(unit_post[[season_var]])
  post_seasons <- post_seasons[!is.na(post_seasons)]
  uncovered <- setdiff(post_seasons, pre_seasons)

  if (length(uncovered) > 0L) {
    stop_lwdid(
      message = sprintf(
        "Unit %s: post-treatment %s(s) %s not in pre-treatment %ss %s.",
        as.character(uid), freq_label,
        paste(sort(uncovered), collapse = ", "),
        freq_label,
        paste(sort(pre_seasons), collapse = ", ")
      ),
      class = c("lwdid_insufficient_quarter_diversity",
                 "lwdid_insufficient_data"),
      unit = uid, missing_quarters = sort(uncovered),
      pre_quarters = sort(pre_seasons),
      post_quarters = sort(post_seasons)
    )
  }
}

# ============================================================================
# Utility: Floating-point safe cohort comparison
# ============================================================================

#' Create boolean mask identifying units belonging to a specific cohort
#'
#' Uses COHORT_FLOAT_TOLERANCE for floating-point comparison.
#' Matches Python get_cohort_mask().
#'
#' @param unit_gvar numeric vector. Unit-level gvar values.
#' @param g numeric. Target cohort (first treatment period).
#' @return logical vector. TRUE for units in cohort g.
#' @keywords internal
.get_cohort_mask <- function(unit_gvar, g) {
  abs(unit_gvar - g) < COHORT_FLOAT_TOLERANCE
}

# ============================================================================
# Utility: Frequency auto-detection
# ============================================================================

#' Detect data frequency from time variable
#'
#' Analyzes time variable to determine frequency (annual/quarterly/monthly/weekly).
#' Uses observations-per-year heuristic for year-like values, interval analysis otherwise.
#' Matches Python detect_frequency() / _detect_frequency_numeric().
#'
#' @param dt data.table
#' @param tvar character(1). Time variable column name.
#' @param ivar character(1) or NULL. Unit identifier for per-unit analysis.
#' @return list with frequency, Q, confidence, method, details
#' @keywords internal
.detect_frequency <- function(dt, tvar, ivar = NULL) {
  result <- list(
    frequency = NULL, Q = NULL, confidence = 0,
    method = NULL, details = list()
  )

  if (!tvar %in% names(dt)) {
    warn_lwdid(
      message = sprintf(
        "Time variable '%s' not found. Cannot detect frequency.", tvar
      ),
      class = "lwdid_data",
      detail = "frequency detection failed",
      action_taken = "skipped frequency detection"
    )
    return(result)
  }

  time_values <- dt[[tvar]]
  time_values <- time_values[!is.na(time_values)]
  if (length(time_values) < 2L) {
    warn_lwdid(
      message = "Insufficient time values for frequency detection (need >= 2).",
      class = "lwdid_data",
      detail = "too few time values",
      action_taken = "skipped frequency detection"
    )
    return(result)
  }

  # Numeric time variable detection
  .detect_frequency_numeric(dt, tvar, ivar)
}

#' Detect frequency from numeric time variable
#'
#' Uses observations-per-year and value range analysis.
#' Matches Python _detect_frequency_numeric().
#'
#' @param dt data.table
#' @param tvar character(1)
#' @param ivar character(1) or NULL
#' @return list with frequency, Q, confidence, method, details
#' @keywords internal
.detect_frequency_numeric <- function(dt, tvar, ivar = NULL) {
  result <- list(
    frequency = NULL, Q = NULL, confidence = 0,
    method = "heuristic", details = list()
  )

  time_values <- dt[[tvar]]
  time_values <- time_values[!is.na(time_values)]
  t_min <- min(time_values)
  t_max <- max(time_values)
  t_range <- t_max - t_min

  result$details$t_min <- t_min
  result$details$t_max <- t_max
  result$details$t_range <- t_range

  looks_like_years <- t_min >= 1900 && t_min <= 2100 &&
    t_max >= 1900 && t_max <= 2100

  if (looks_like_years && t_range > 0) {
    # Count observations per year per unit
    if (!is.null(ivar) && ivar %in% names(dt)) {
      obs_per_year_list <- numeric(0)
      n_units <- 0L
      for (uid in unique(dt[[ivar]])) {
        unit_times <- dt[dt[[ivar]] == uid][[tvar]]
        unit_times <- unit_times[!is.na(unit_times)]
        if (length(unit_times) >= 2L) {
          year_counts <- table(as.integer(unit_times))
          obs_per_year_list <- c(obs_per_year_list,
                                 stats::median(year_counts))
          n_units <- n_units + 1L
        }
      }
      result$details$n_units_analyzed <- n_units
      obs_per_year <- if (length(obs_per_year_list) > 0)
        stats::median(obs_per_year_list) else NULL
    } else {
      year_counts <- table(as.integer(time_values))
      obs_per_year <- if (length(year_counts) > 0)
        stats::median(year_counts) else NULL
      result$details$n_units_analyzed <- 1L
    }

    if (!is.null(obs_per_year)) {
      result$details$obs_per_year <- obs_per_year
      result$method <- "obs_per_year"

      # Frequency mapping: (expected, tolerance, name, Q)
      freq_mapping <- list(
        list(1, 0.5, "annual", 1L),
        list(4, 1.5, "quarterly", 4L),
        list(12, 3, "monthly", 12L),
        list(52, 10, "weekly", 52L)
      )

      best_distance <- Inf
      for (fm in freq_mapping) {
        distance <- abs(obs_per_year - fm[[1]])
        if (distance < fm[[2]] && distance < best_distance) {
          best_distance <- distance
          result$frequency <- fm[[3]]
          result$Q <- fm[[4]]
          result$confidence <- min(1, 1 - distance / fm[[2]])
        }
      }
    }
  } else {
    # Non-year numeric: interval analysis
    if (!is.null(ivar) && ivar %in% names(dt)) {
      intervals <- numeric(0)
      for (uid in unique(dt[[ivar]])) {
        unit_times <- sort(dt[dt[[ivar]] == uid][[tvar]])
        unit_times <- unit_times[!is.na(unit_times)]
        if (length(unit_times) >= 2L) {
          intervals <- c(intervals, diff(unit_times))
        }
      }
    } else {
      sorted_times <- sort(time_values)
      intervals <- diff(sorted_times)
    }

    if (length(intervals) > 0L) {
      median_interval <- stats::median(intervals)
      result$details$median_interval <- median_interval
      result$method <- "interval"

      if (median_interval == 1) {
        result$confidence <- 0.3
        warn_lwdid(
          message = paste0(
            "Time variable appears to be consecutive integer index. ",
            "Cannot reliably detect frequency. ",
            "Please specify Q explicitly."
          ),
          class = "lwdid_data",
          detail = "ambiguous frequency",
          action_taken = "frequency detection inconclusive"
        )
      }
    }
  }

  result
}

# ============================================================================
# Deferred cross-parameter checks (Layer 10)
# ============================================================================

#' Execute deferred cross-parameter checks requiring stag_info
#'
#' These checks are logically part of Layer 6 but deferred to Layer 10
#' because they require cohort/NT information from Layer 9.
#'
#' @param stag_info list. Output from .validate_staggered_data().
#' @param aggregate character(1). Aggregation type.
#' @param control_group character(1). Control group strategy.
#' @param gid scalar or NULL. Unit ID to highlight.
#' @param dt data.table. Cleaned data.
#' @param ivar character(1). Unit identifier column name.
#' @param gvar character(1). Cohort variable column name.
#' @return character(1). Resolved control_group.
#' @keywords internal
.validate_deferred_staggered_checks <- function(
  stag_info, aggregate, control_group, gid, dt, ivar, gvar
) {
  has_nt <- stag_info$has_never_treated
  cg_resolved <- control_group

  # J7: aggregate requires NT units
  if (aggregate %in% c("cohort", "overall") && !has_nt) {
    stop_lwdid(
      message = paste0(
        "aggregate='", aggregate, "' requires never-treated (NT) units for ",
        "cross-cohort aggregation, but no NT units found in data.\n",
        "Use aggregate='none' to estimate (g,r)-specific effects only."
      ),
      class = "lwdid_no_never_treated",
      n_never_treated = 0L,
      detail = paste0("aggregate='", aggregate, "' requires NT")
    )
  }

  # J8: never_treated control group requires NT units
  if (control_group == "never_treated" && !has_nt) {
    stop_lwdid(
      message = paste0(
        "control_group='never_treated' specified but no never-treated ",
        "units found.\n",
        "Use control_group='not_yet_treated' or control_group='auto' instead."
      ),
      class = "lwdid_no_never_treated",
      n_never_treated = 0L,
      detail = "control_group='never_treated' requires NT"
    )
  }

  # J9: aggregate + non-NT control_group → auto-switch to never_treated

  # Both not_yet_treated and all_others need switching for cohort/overall
  # aggregation (FATAL-004 defense). all_others includes already-treated

  # units which violates the common control group requirement.
  if (aggregate %in% c("cohort", "overall") &&
      control_group %in% c("not_yet_treated", "all_others") && has_nt) {
    warn_lwdid(
      message = paste0(
        "aggregate='", aggregate,
        "' with control_group='", control_group, "': ",
        "auto-switching to control_group='never_treated' for aggregation. ",
        "Cross-cohort aggregation requires a common control group (NT units)."
      ),
      class = "lwdid_control_group_switch",
      detail = "auto-switch control_group for aggregation",
      action_taken = "control_group switched to 'never_treated'"
    )
    cg_resolved <- "never_treated"
  }

  # J10: all_others bias warning for non-aggregation scenarios
  if (control_group == "all_others" &&
      aggregate %in% c("none", "event_time")) {
    warn_lwdid(
      paste("control_group = 'all_others' includes already-treated units",
            "as controls, which may introduce bias. Consider using",
            "'not_yet_treated' instead."),
      class = "lwdid_data",
      detail = "all_others_bias"
    )
  }

  # control_group = "auto" resolution
  if (control_group == "auto") {
    cg_resolved <- if (has_nt) "never_treated" else "not_yet_treated"
  }

  # gid validation: must be a treated (non-NT) unit
  if (!is.null(gid)) {
    unit_gvar <- dt[, .(g = .SD[[1L]][1L]), by = c(ivar), .SDcols = gvar]
    gid_row <- unit_gvar[unit_gvar[[ivar]] == gid]
    if (nrow(gid_row) == 0L) {
      stop_lwdid(
        message = sprintf(
          "gid=%s not found in data. Available unit IDs: %s",
          as.character(gid),
          paste(head(unique(dt[[ivar]]), 10), collapse = ", ")
        ),
        class = "lwdid_invalid_parameter",
        param = "gid", value = gid, allowed = "existing unit ID"
      )
    }
    if (is_never_treated(gid_row$g[1])) {
      stop_lwdid(
        message = sprintf(
          "gid=%s is a never-treated unit. gid must be a treated unit.",
          as.character(gid)
        ),
        class = "lwdid_invalid_parameter",
        param = "gid", value = gid, allowed = "treated unit ID"
      )
    }
  }

  cg_resolved
}

#' Validate gid refers to a treated unit (Common Timing mode)
#'
#' @param gid scalar or NULL. Unit ID to validate.
#' @param dt data.table with d_ column.
#' @param ivar character(1). Unit identifier column name.
#' @param treat_col character(1). Treatment indicator column ("d_").
#' @return NULL (invisible).
#' @keywords internal
.validate_gid_treated <- function(gid, dt, ivar, treat_col) {
  if (is.null(gid)) return(invisible(NULL))

  if (!gid %in% dt[[ivar]]) {
    stop_lwdid(
      message = sprintf(
        "gid=%s not found in data. Available unit IDs: %s",
        as.character(gid),
        paste(head(unique(dt[[ivar]]), 10), collapse = ", ")
      ),
      class = "lwdid_invalid_parameter",
      param = "gid", value = gid, allowed = "existing unit ID"
    )
  }

  unit_d <- dt[dt[[ivar]] == gid,
               .(d_max = max(.SD[[1L]], na.rm = TRUE)), .SDcols = treat_col]
  if (nrow(unit_d) == 0L || unit_d$d_max != 1L) {
    stop_lwdid(
      message = sprintf(
        "gid=%s is not a treated unit (d==0). gid must be a treated unit (d==1).",
        as.character(gid)
      ),
      class = "lwdid_invalid_parameter",
      param = "gid", value = gid,
      allowed = "treated unit ID (d==1)"
    )
  }
  invisible(NULL)
}


# ============================================================================
# Main entry point: validate_inputs()
# ============================================================================

#' Validate and prepare inputs for lwdid estimation
#'
#' Single entry point for all input validation. Performs a 10-layer
#' validation pipeline and returns a standardized result list.
#'
#' @param data data.frame or data.table. Panel data in long format.
#' @param y character(1). Outcome variable column name.
#' @param d character(1) or NULL. Treatment indicator (Common Timing).
#' @param ivar character(1). Unit identifier column name.
#' @param tvar character(1) or character(2). Time variable(s).
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param rolling character(1). Transformation method.
#' @param gvar character(1) or NULL. Cohort variable (Staggered).
#' @param control_group character(1). Control group strategy.
#' @param aggregate character(1). Aggregation type.
#' @param estimator character(1). Estimator type.
#' @param controls character vector or NULL. Control variables.
#' @param ps_controls character vector or NULL. PS model controls.
#' @param trim_threshold numeric(1). PS trimming threshold.
#' @param n_neighbors integer(1). PSM neighbors.
#' @param caliper numeric(1) or NULL. PSM caliper.
#' @param with_replacement logical(1). PSM with replacement.
#' @param match_order character(1). PSM match order.
#' @param return_diagnostics logical(1). Return PS diagnostics.
#' @param vce character(1) or NULL. Variance estimator type.
#' @param cluster_var character(1) or NULL. Cluster variable.
#' @param alpha numeric(1). Significance level.
#' @param ri logical(1). Randomization inference.
#' @param rireps integer(1). RI replications.
#' @param seed integer(1) or NULL. Random seed.
#' @param ri_method character(1). RI method.
#' @param season_var character(1) or NULL. Seasonal indicator.
#' @param Q integer(1). Seasons per cycle.
#' @param auto_detect_frequency logical(1). Auto-detect frequency.
#' @param include_pretreatment logical(1). Include pretreatment.
#' @param pretreatment_test logical(1). Pretreatment test.
#' @param pretreatment_alpha numeric(1). Pretreatment alpha.
#' @param exclude_pre_periods integer(1). Periods to exclude.
#' @param balanced_panel character(1). Balance handling.
#' @param graph logical(1). Generate graph.
#' @param gid NULL or scalar. Unit ID to highlight.
#' @param graph_options list or NULL. Graph options.
#' @param verbose character(1). Verbosity level.
#' @param att_gt logical(1). ATT(g,t) output.
#' @return A named list with validated data and metadata.
#' @keywords internal

validate_inputs <- function(
  data,
  y,
  d = NULL,
  ivar,
  tvar,
  post = NULL,
  rolling = "demean",
  gvar = NULL,
  control_group = "not_yet_treated",
  aggregate = "cohort",
  estimator = "ra",
  controls = NULL,
  ps_controls = NULL,
  trim_threshold = 0.01,
  trim_method = "clip",
  n_neighbors = 1L,
  caliper = NULL,
  caliper_scale = "sd",
  with_replacement = TRUE,
  match_order = "data",
  se_method = NULL,
  boot_reps = 200L,
  return_diagnostics = FALSE,
  vce = NULL,
  cluster_var = NULL,
  alpha = 0.05,
  ri = FALSE,
  rireps = 1000L,
  seed = NULL,
  ri_method = "bootstrap",
  season_var = NULL,
  Q = 4L,
  auto_detect_frequency = FALSE,
  include_pretreatment = FALSE,
  pretreatment_test = TRUE,
  pretreatment_alpha = 0.05,
  exclude_pre_periods = 0L,
  balanced_panel = "warn",
  graph = FALSE,
  gid = NULL,
  graph_options = NULL,
  verbose = "default",
  att_gt = FALSE,
  registry = NULL
) {

  # ==== Layer 1: Reserved columns + empty data ====
  .validate_reserved_columns(data)
  .validate_data_not_empty(data)

  # ==== Layer 2: Parameter type/range validation ====
  .validate_string_param(y, "y", data)
  .validate_string_param(ivar, "ivar", data)
  .validate_tvar(tvar, data)
  .validate_optional_string_param(d, "d", data)
  .validate_optional_string_param(post, "post", data)
  .validate_optional_string_param(gvar, "gvar", data)
  rolling_lower <- .validate_rolling_parameter(
    rolling, tvar, season_var
  )
  .validate_choice(estimator, VALID_ESTIMATORS, "estimator")
  .validate_choice(control_group, VALID_CONTROL_GROUPS,
                   "control_group")
  .validate_choice(aggregate, VALID_AGGREGATES, "aggregate")
  .validate_choice(balanced_panel, VALID_BALANCED_PANEL,
                   "balanced_panel")
  .validate_choice(verbose, VALID_VERBOSE, "verbose")
  .validate_choice(match_order, VALID_MATCH_ORDER, "match_order")
  .validate_choice(trim_method, VALID_TRIM_METHODS, "trim_method")
  .validate_choice(caliper_scale, VALID_CALIPER_SCALES,
                   "caliper_scale")
  .validate_choice(ri_method, VALID_RI_METHODS, "ri_method")
  .validate_vce(vce)

  .validate_numeric_range(alpha, 0, 1, "alpha", exclusive = TRUE)
  .validate_numeric_range(pretreatment_alpha, 0, 1,
                          "pretreatment_alpha", exclusive = TRUE)
  .validate_numeric_range(trim_threshold, 0, 0.5,
                          "trim_threshold", exclusive = TRUE)
  .validate_positive_integer(n_neighbors, "n_neighbors")
  .validate_positive_integer(Q, "Q", min_val = 2L)
  .validate_positive_integer(rireps, "rireps")
  .validate_positive_integer(boot_reps, "boot_reps")
  .validate_nonneg_integer(exclude_pre_periods,
                           "exclude_pre_periods")
  .validate_logical(ri, "ri")
  .validate_logical(graph, "graph")
  .validate_logical(att_gt, "att_gt")
  .validate_logical(with_replacement, "with_replacement")
  .validate_logical(include_pretreatment, "include_pretreatment")
  .validate_logical(pretreatment_test, "pretreatment_test")
  .validate_logical(auto_detect_frequency,
                    "auto_detect_frequency")
  .validate_logical(return_diagnostics, "return_diagnostics")
  .validate_optional_positive_numeric(caliper, "caliper")
  if (!is.null(se_method)) {
    if (!is.character(se_method) || length(se_method) != 1L) {
      stop_lwdid(
        sprintf("se_method must be a single character string or NULL, got: %s",
                deparse(se_method)),
        class = "lwdid_invalid_parameter",
        param = "se_method", value = se_method
      )
    }
  }
  .validate_optional_integer(seed, "seed")
  .validate_controls_columns(controls, data)
  .validate_controls_columns(ps_controls, data)
  .validate_optional_string_param(season_var, "season_var", data)
  .validate_season_var_values(data, season_var, Q)
  .validate_optional_string_param(cluster_var, "cluster_var", data)
  .validate_graph_options(graph_options)

  # ==== Layer 3a: Mode identification (no data modification) ====
  mode <- .identify_mode(d, post, gvar)

  # ==== Layer 4: Data type validation ====
  .validate_outcome_dtype(data, y)
  .validate_controls_dtype(data, controls)

  # ==== Layer 5: Time-invariance (CT only) ====
  if (mode == "common_timing") {
    .validate_treatment_time_invariance(data, d, ivar)
    .validate_time_invariant_controls(data, ivar, controls)
  }

  # ==== Layer 6: Cross-parameter consistency ====
  .validate_cross_param_consistency(
    vce = vce, cluster_var = cluster_var,
    rolling = rolling_lower, season_var = season_var,
    tvar = tvar, estimator = estimator,
    controls = controls, ps_controls = ps_controls,
    aggregate = aggregate, control_group = control_group,
    graph = graph, ri = ri, rireps = rireps,
    ri_method = ri_method,
    exclude_pre_periods = exclude_pre_periods,
    Q = Q, mode = mode
  )

  # ==== Layer 7: Data format conversion (copy semantics) ====
  dt <- .convert_to_datatable(data)
  id_result <- .convert_string_id(dt, ivar)
  dt <- id_result$data
  id_mapping <- id_result$id_mapping

  # ==== Layer 3b: Binarization on data.table copy ====
  if (mode == "common_timing") {
    d_vec <- dt[[d]]
    post_vec <- dt[[post]]
    dt[, d_ := .binarize_with_na(d_vec)]
    dt[, post_ := .binarize_with_na(post_vec)]
  }

  # ==== Layer 8: Missing value handling ====
  dt <- .handle_missing_values(
    dt, y, ivar, tvar, mode, d, post
  )

  # ==== Layer 9: Time index creation ====
  ti_result <- .create_time_index(dt, tvar)
  dt <- ti_result$data
  is_quarterly <- ti_result$is_quarterly

  # ==== Layer 9 (cont): Staggered data validation ====
  stag_info <- NULL
  if (mode == "staggered") {
    stag_info <- .validate_staggered_data(
      dt, gvar, ivar, tvar, y
    )
    # Emit collected staggered warnings
    for (w in stag_info$warning_list) {
      warn_lwdid(
        message = w,
        class = "lwdid_data",
        detail = "staggered data diagnostic",
        action_taken = "continued"
      )
    }
  }

  # ==== Layer 10: Panel structure validation ====
  # --- Duplicate check ---
  .validate_no_duplicates(dt, ivar, tvar)

  # --- Time continuity ---
  T_total <- as.integer(max(dt$tindex))
  unique_tindex <- sort(unique(dt$tindex))
  expected_tindex <- seq_len(T_total)
  if (!identical(as.integer(unique_tindex), expected_tindex)) {
    missing_periods <- setdiff(expected_tindex, unique_tindex)
    stop_lwdid(
      message = sprintf(
        "Time series is discontinuous. Missing periods: %s",
        paste(missing_periods, collapse = ", ")
      ),
      class = "lwdid_time_discontinuity",
      gaps = missing_periods, ivar = ivar,
      tvar = if (length(tvar) == 1L) tvar else tvar[1]
    )
  }

  # --- Common Timing specific checks ---
  K <- NULL
  tpost1 <- NULL
  n_treated <- NULL
  n_control <- NULL
  n_pre <- NULL

  if (mode == "common_timing") {
    ct_info <- .validate_common_timing_structure(dt, ivar)
    K <- ct_info$K
    tpost1 <- ct_info$tpost1
    n_treated <- ct_info$n_treated
    n_control <- ct_info$n_control
    n_pre <- K  # S-1 = K = tpost1 - 1
  }

  # --- Sample size ---
  n_units <- data.table::uniqueN(dt[[ivar]])
  if (n_units < 3L) {
    stop_lwdid(
      message = sprintf(
        "Insufficient sample size: N=%d (need N >= 3).",
        n_units
      ),
      class = "lwdid_insufficient_data",
      n = n_units, n_treated = 0L, n_control = 0L
    )
  }

  # --- Small sample warnings ---
  if (n_units <= 30L) {
    warn_lwdid(
      message = sprintf(
        "Small sample (N=%d). Consider checking normality assumptions.",
        n_units
      ),
      class = "lwdid_small_sample",
      n = n_units, n_treated = 0L, n_control = 0L
    )
  }

  # --- Panel balance ---
  balanced <- .validate_panel_balance(
    dt, ivar, balanced_panel, mode, gvar, tvar
  )

  # --- Deferred cross-parameter checks (Layer 10) ---
  cg_resolved <- control_group
  if (mode == "staggered" && !is.null(stag_info)) {
    cg_resolved <- .validate_deferred_staggered_checks(
      stag_info, aggregate, control_group, gid, dt, ivar, gvar
    )
    n_treated <- stag_info$n_treated
    n_control <- if (stag_info$has_never_treated) {
      stag_info$n_never_treated
    } else {
      n_units - stag_info$n_treated
    }
  }
  if (mode == "common_timing") {
    .validate_gid_treated(gid, dt, ivar, "d_")
    # Resolve auto control_group for CT (not applicable, but safe)
    if (control_group == "auto") {
      cg_resolved <- "not_yet_treated"
    }
  }

  # --- Pre-period sufficiency ---
  cohorts_for_check <- if (mode == "staggered" &&
                           !is.null(stag_info)) {
    stag_info$cohorts
  } else {
    NULL
  }
  T_min_val <- if (mode == "staggered" && !is.null(stag_info)) {
    stag_info$T_min
  } else {
    tvar_col <- if (length(tvar) == 1L) tvar else tvar[1]
    as.integer(min(dt[[tvar_col]], na.rm = TRUE))
  }
  season_col <- season_var
  if (is.null(season_col) && length(tvar) == 2L) {
    season_col <- tvar[2]
  }

  seasonal_gate_applied <- FALSE
  if (mode == "common_timing" &&
      rolling_lower %in% c("demeanq", "detrendq") &&
      !is.null(season_col)) {
    .validate_seasonal_inputs(
      data = dt,
      y = y,
      season_var = season_col,
      Q = Q,
      ivar = ivar,
      tindex = "tindex",
      post = "post_",
      method = rolling_lower,
      exclude_pre_periods = exclude_pre_periods,
      min_global_pre_periods = if (rolling_lower == "detrendq") 2L else 1L,
      registry = registry
    )
    seasonal_gate_applied <- TRUE
  } else {
    .validate_pre_period_sufficiency(
      mode = mode, rolling = rolling_lower,
      n_pre_periods = if (!is.null(n_pre)) n_pre else 0L,
      exclude_pre_periods = exclude_pre_periods,
      Q = Q, cohorts = cohorts_for_check,
      T_min = T_min_val, dt = dt, ivar = ivar, gvar = gvar
    )
  }

  # --- Seasonal diversity (if demeanq/detrendq) ---
  if (rolling_lower %in% c("demeanq", "detrendq") &&
      !seasonal_gate_applied &&
      !is.null(season_col)) {
    if (mode == "common_timing") {
      .validate_season_diversity(
        dt, ivar, season_col, "post_", Q
      )
    } else {
      .validate_season_diversity(
        dt, ivar, season_col, NULL, Q,
        gvar = gvar, tvar = tvar
      )
    }
  }

  # --- Frequency detection ---
  freq_info <- NULL
  if (isTRUE(auto_detect_frequency)) {
    tvar_col <- if (length(tvar) == 1L) tvar else tvar[1]
    freq_info <- .detect_frequency(dt, tvar_col, ivar)
  }

  if (rolling_lower %in% c("demeanq", "detrendq") && !is.null(season_col)) {
    attr(dt, "lwdid_seasonal_config") <- list(
      post = "post_",
      season_var = season_col,
      Q = as.integer(Q)
    )
  }

  # ==== Build return value ====
  n_obs <- nrow(dt)
  n_periods <- data.table::uniqueN(dt$tindex)
  periods <- sort(unique(dt$tindex))

  # Time range
  tvar_col <- if (length(tvar) == 1L) tvar else tvar[1]
  T_min_out <- as.integer(min(dt[[tvar_col]], na.rm = TRUE))
  T_max_out <- as.integer(max(dt[[tvar_col]], na.rm = TRUE))

  # Staggered-specific fields
  cohorts_out <- NULL
  n_cohorts_out <- NULL
  cohort_sizes_out <- NULL
  has_nt <- NA
  n_nt <- 0L
  cohort_pre <- NULL

  if (mode == "staggered" && !is.null(stag_info)) {
    cohorts_out <- stag_info$cohorts
    n_cohorts_out <- stag_info$n_cohorts
    cohort_sizes_out <- stag_info$cohort_sizes
    has_nt <- stag_info$has_never_treated
    n_nt <- stag_info$n_never_treated
    # Per-cohort pre-period counts
    cohort_pre <- setNames(
      as.list(cohorts_out - T_min_out),
      as.character(cohorts_out)
    )
  }

  list(
    mode = mode,
    data = dt,
    n_obs = n_obs,
    n_units = n_units,
    n_periods = n_periods,
    K = K,
    tpost1 = tpost1,
    T_min = T_min_out,
    T_max = T_max_out,
    is_quarterly = is_quarterly,
    n_pre_periods = n_pre,
    cohort_pre_periods = cohort_pre,
    cohorts = cohorts_out,
    n_cohorts = n_cohorts_out,
    cohort_sizes = cohort_sizes_out,
    has_nt = has_nt,
    n_treated = n_treated,
    n_control = n_control,
    n_nt = n_nt,
    periods = periods,
    balanced = balanced,
    id_mapping = id_mapping,
    control_group_resolved = cg_resolved,
    freq_info = freq_info,

    validated_params = list(
      depvar = y,
      ivar = ivar,
      tvar = tvar,
      d = d,
      post = post,
      gvar = gvar,
      rolling = rolling_lower,
      estimator = estimator,
      controls = controls,
      ps_controls = ps_controls,
      trim_threshold = trim_threshold,
      trim_method = trim_method,
      n_neighbors = n_neighbors,
      caliper = caliper,
      caliper_scale = caliper_scale,
      with_replacement = with_replacement,
      match_order = match_order,
      se_method = se_method,
      boot_reps = boot_reps,
      return_diagnostics = return_diagnostics,
      vce = vce,
      cluster_var = cluster_var,
      alpha = alpha,
      ri = ri,
      rireps = rireps,
      seed = seed,
      ri_method = ri_method,
      season_var = season_var,
      Q = Q,
      auto_detect_frequency = auto_detect_frequency,
      include_pretreatment = include_pretreatment,
      pretreatment_test = pretreatment_test,
      pretreatment_alpha = pretreatment_alpha,
      exclude_pre_periods = exclude_pre_periods,
      balanced_panel = balanced_panel,
      graph = graph,
      gid = gid,
      graph_options = graph_options,
      verbose = verbose,
      att_gt = att_gt,
      aggregate = aggregate,
      control_group = control_group
    )
  )
}
