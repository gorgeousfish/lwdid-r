#' @title lwdid Custom Condition System
#' @description Structured error and warning conditions for the lwdid package.
#'   Implements 18 error conditions and 5 warning conditions mirroring the
#'   Python lwdid-py v0.2.3 exception hierarchy.
#' @param param Character scalar naming the parameter that failed validation.
#' @param value Observed invalid value.
#' @param allowed Allowed values or a human-readable constraint summary.
#' @param call Call object associated with the condition, or NULL.
#' @param n,n_treated,n_control Integer counts for total, treated, and control
#'   sample sizes.
#' @param min_obs,max_obs Integer minimum and maximum observations per unit.
#' @param n_incomplete_units Integer number of incomplete units.
#' @param pct_unbalanced Numeric percentage of incomplete units.
#' @param gaps Vector of detected time gaps.
#' @param ivar,tvar Character names of the unit and time variables.
#' @param column Character name of a missing required column.
#' @param available Character vector of available columns.
#' @param reps_completed Integer number of completed replications before
#'   failure.
#' @param error_detail Character summary of the triggering error.
#' @param plot_type Character plot type associated with a visualization error.
#' @param detail Character detail describing the issue.
#' @param gvar Character name of the cohort variable.
#' @param invalid_values Vector of invalid cohort values.
#' @param aggregate_type Character aggregation type associated with the error.
#' @param method Character transformation or rolling method.
#' @param vce_type Character variance estimator type.
#' @param treat_var Character treatment-variable name.
#' @param n_units Integer number of units.
#' @param required Integer minimum required count.
#' @param rolling Character rolling method.
#' @param cohort Cohort identifier associated with the condition, if any.
#' @param excluded Integer number of excluded periods, if any.
#' @param unit Unit identifier associated with the condition.
#' @param missing_quarters Integer vector of missing quarters.
#' @param pre_quarters,post_quarters Integer vectors of observed pre/post
#'   quarters.
#' @param aggregate Character aggregation choice.
#' @param control_group Character control-group choice.
#' @param violation Character summary of the violated aggregation rule.
#' @param min_cell_size Integer minimum acceptable cell size.
#' @param max_observed_size Integer largest observed cell size.
#' @param n_cells Integer number of cells evaluated.
#' @param period Period identifier associated with the warning, if any.
#' @param n_extreme Integer count of extreme propensity scores.
#' @param ps_range Numeric length-2 vector giving the propensity-score range.
#' @param max_weight Numeric maximum inverse-probability weight.
#' @param source_function Character function name that emitted the numerical
#'   warning.
#' @param condition_number Numeric condition number, if available.
#' @param action_taken Character remediation or fallback that was applied.
#' @param model Character model name associated with a convergence warning.
#' @param iterations,max_iterations Integer iteration counts.
#' @param tolerance Numeric convergence tolerance.
#' @param original,switched_to Character original and replacement choices in an
#'   automatic control-group switch.
#' @name lwdid-conditions
#' @family lwdid-conditions
NULL

# ============================================================================
# Core Factory Functions
# ============================================================================

#' Create an lwdid error condition object
#'
#' @param message Character string describing the error.
#' @param class Character vector of S3 class names (specific to general).
#'   `"lwdid_error"`, `"error"`, `"condition"` are appended automatically.
#' @param ... Additional named fields attached to the condition object.
#' @param call The call associated with the condition, or NULL.
#' @return A condition object (not thrown).
#' @export
#' @family lwdid-conditions
#' @examples
#' cond <- lwdid_error("something went wrong", class = "lwdid_missing_column",
#'                     column = "y")
#' class(cond)
lwdid_error <- function(message, class, ..., call = NULL) {
  if ("lwdid_invalid_param" %in% class &&
      !"lwdid_invalid_parameter" %in% class) {
    class <- c("lwdid_invalid_param", "lwdid_invalid_parameter",
               setdiff(class, "lwdid_invalid_param"))
  }
  tail_classes <- c("lwdid_error", "error", "condition")
  full_class <- unique(c(class, tail_classes))
  structure(
    class = full_class,
    list(message = message, call = call, ...)
  )
}

#' Create an lwdid warning condition object
#'
#' @param message Character string describing the warning.
#' @param class Character vector of S3 class names (specific to general).
#'   `"lwdid_warning"`, `"warning"`, `"condition"` are appended automatically.
#' @param ... Additional named fields attached to the condition object.
#' @param call The call associated with the condition, or NULL.
#' @return A condition object (not thrown).
#' @export
#' @family lwdid-conditions
lwdid_warning <- function(message, class, ..., call = NULL) {
  tail_classes <- c("lwdid_warning", "warning", "condition")
  full_class <- unique(c(class, tail_classes))
  structure(
    class = full_class,
    list(message = message, call = call, ...)
  )
}

# ============================================================================
# Convenience Throw Functions
# ============================================================================

#' Construct and immediately throw an lwdid error
#'
#' @inheritParams lwdid_error
#' @return Does not return; throws an error condition.
#' @export
#' @family lwdid-conditions
stop_lwdid <- function(message, class, ..., call = sys.call(-1)) {
  cond <- lwdid_error(message, class, ..., call = call)
  stop(cond)
}

#' Construct and immediately throw an lwdid warning
#'
#' @inheritParams lwdid_warning
#' @return Returns NULL invisibly after issuing the warning.
#' @export
#' @family lwdid-conditions
warn_lwdid <- function(message, class, ..., call = sys.call(-1)) {
  cond <- lwdid_warning(message, class, ..., call = call)
  warning(cond)
  invisible(NULL)
}

# ============================================================================
# Error Condition Convenience Constructors (18 types)
# ============================================================================

# --- No intermediate parent (direct -> lwdid_error) ---

#' @rdname lwdid-conditions
#' @export
lwdid_invalid_parameter_error <- function(param, value, allowed, call = NULL) {
  msg <- sprintf(
    "Invalid parameter '%s': got '%s'. Allowed values: %s",
    param, as.character(value), paste(allowed, collapse = ", ")
  )
  lwdid_error(msg, class = "lwdid_invalid_parameter",
              param = param, value = value, allowed = allowed, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_insufficient_data_error <- function(n, n_treated, n_control, call = NULL) {
  msg <- sprintf(
    "Insufficient data for estimation: N=%d (treated=%d, control=%d). Minimum required: 3.",
    n, n_treated, n_control
  )
  lwdid_error(msg, class = "lwdid_insufficient_data",
              n = n, n_treated = n_treated, n_control = n_control, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_unbalanced_panel_error <- function(min_obs, max_obs,
                                          n_incomplete_units, pct_unbalanced,
                                          call = NULL) {
  msg <- sprintf(
    "Unbalanced panel detected: obs per unit ranges from %d to %d. %d units (%.1f%%) are incomplete.",
    min_obs, max_obs, n_incomplete_units, pct_unbalanced
  )
  lwdid_error(msg, class = "lwdid_unbalanced_panel",
              min_obs = min_obs, max_obs = max_obs,
              n_incomplete_units = n_incomplete_units,
              pct_unbalanced = pct_unbalanced, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_time_discontinuity_error <- function(gaps, ivar, tvar, call = NULL) {
  msg <- sprintf(
    "Time discontinuity detected in '%s' (unit var: '%s'): gaps at %s.",
    tvar, ivar, paste(gaps, collapse = ", ")
  )
  lwdid_error(msg, class = "lwdid_time_discontinuity",
              gaps = gaps, ivar = ivar, tvar = tvar, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_missing_column_error <- function(column, available, call = NULL) {
  msg <- sprintf(
    "Required column '%s' not found in data. Available columns: %s",
    column, paste(available, collapse = ", ")
  )
  lwdid_error(msg, class = "lwdid_missing_column",
              column = column, available = available, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_randomization_error_cond <- function(reps_completed, error_detail,
                                            call = NULL) {
  msg <- sprintf(
    "Randomization inference failed after %d replications: %s",
    reps_completed, error_detail
  )
  lwdid_error(msg, class = "lwdid_randomization_error",
              reps_completed = reps_completed, error_detail = error_detail,
              call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_visualization_error_cond <- function(plot_type, detail, call = NULL) {
  msg <- sprintf("Visualization error for '%s' plot: %s", plot_type, detail)
  lwdid_error(msg, class = "lwdid_visualization_error",
              plot_type = plot_type, detail = detail, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_invalid_staggered_data_error <- function(gvar, invalid_values, detail,
                                                call = NULL) {
  msg <- sprintf(
    "Invalid staggered data in '%s': %s. Invalid values: %s",
    gvar, detail, paste(head(invalid_values, 10), collapse = ", ")
  )
  lwdid_error(msg, class = "lwdid_invalid_staggered_data",
              gvar = gvar, invalid_values = invalid_values,
              detail = detail, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_aggregation_error_cond <- function(aggregate_type, detail, call = NULL) {
  msg <- sprintf("Aggregation failed for type '%s': %s", aggregate_type, detail)
  lwdid_error(msg, class = "lwdid_aggregation_error",
              aggregate_type = aggregate_type, detail = detail, call = call)
}

# --- With intermediate parent lwdid_invalid_parameter ---

#' @rdname lwdid-conditions
#' @export
lwdid_invalid_rolling_error <- function(method, allowed = c("demean", "detrend", "demeanq", "detrendq"), call = NULL) {
  msg <- sprintf(
    "Invalid rolling method '%s'. Must be one of: %s",
    method, paste(allowed, collapse = ", ")
  )
  lwdid_error(msg,
              class = c("lwdid_invalid_rolling", "lwdid_invalid_parameter"),
              method = method, allowed = allowed, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_invalid_vce_error <- function(vce_type, allowed = c("NULL", "robust", "hc0", "hc1", "hc2", "hc3", "hc4", "cluster", "bootstrap"), call = NULL) {
  msg <- sprintf(
    "Invalid VCE type '%s'. Must be one of: %s",
    vce_type, paste(allowed, collapse = ", ")
  )
  lwdid_error(msg,
              class = c("lwdid_invalid_vce", "lwdid_invalid_parameter"),
              vce_type = vce_type, allowed = allowed, call = call)
}

# --- With intermediate parent lwdid_insufficient_data ---

#' @rdname lwdid-conditions
#' @export
lwdid_no_treated_error <- function(treat_var, n_units, call = NULL) {
  msg <- sprintf(
    "No treated units found: all %d units have %s = 0.",
    n_units, treat_var
  )
  lwdid_error(msg,
              class = c("lwdid_no_treated", "lwdid_insufficient_data"),
              treat_var = treat_var, n_units = n_units, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_no_control_error <- function(treat_var, n_units, call = NULL) {
  msg <- sprintf(
    "No control units found: all %d units have %s = 1.",
    n_units, treat_var
  )
  lwdid_error(msg,
              class = c("lwdid_no_control", "lwdid_insufficient_data"),
              treat_var = treat_var, n_units = n_units, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_insufficient_pre_periods_error <- function(available, required, rolling,
                                                  cohort = NULL,
                                                  excluded = NULL,
                                                  call = NULL) {
  msg <- sprintf(
    "Insufficient pre-treatment periods for '%s': available=%d, required=%d.",
    rolling, available, required
  )
  if (!is.null(cohort)) {
    msg <- paste0(msg, sprintf(" Cohort %d", cohort))
  }
  if (!is.null(excluded) && excluded > 0L) {
    msg <- paste0(msg, sprintf(", excluded %d periods.", excluded))
  }
  lwdid_error(msg,
              class = c("lwdid_insufficient_pre_periods",
                        "lwdid_insufficient_data"),
              available = available, required = required, rolling = rolling,
              cohort = cohort, excluded = excluded, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_insufficient_quarter_diversity_error <- function(unit, missing_quarters,
                                                        pre_quarters,
                                                        post_quarters,
                                                        call = NULL) {
  msg <- sprintf(
    "Unit %s: post-treatment quarters {%s} not found in pre-treatment quarters {%s}.",
    as.character(unit),
    paste(missing_quarters, collapse = ", "),
    paste(pre_quarters, collapse = ", ")
  )
  lwdid_error(msg,
              class = c("lwdid_insufficient_quarter_diversity",
                        "lwdid_insufficient_data"),
              unit = unit, missing_quarters = missing_quarters,
              pre_quarters = pre_quarters, post_quarters = post_quarters,
              call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_no_never_treated_error <- function(aggregate, control_group,
                                          call = NULL) {
  msg <- sprintf(
    "No never-treated units found. Required for aggregate='%s' with control_group='%s'.",
    aggregate, control_group
  )
  lwdid_error(msg,
              class = c("lwdid_no_never_treated", "lwdid_insufficient_data"),
              aggregate = aggregate, control_group = control_group,
              call = call)
}

# --- With intermediate parent lwdid_aggregation_error ---

#' @rdname lwdid-conditions
#' @export
lwdid_invalid_aggregation_error <- function(aggregate_type, violation, detail,
                                             call = NULL) {
  msg <- sprintf(
    "Invalid aggregation (%s): %s. %s",
    aggregate_type, violation, detail
  )
  lwdid_error(msg,
              class = c("lwdid_invalid_aggregation",
                        "lwdid_aggregation_error"),
              aggregate_type = aggregate_type, violation = violation,
              detail = detail, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_insufficient_cell_size_error <- function(min_cell_size,
                                                max_observed_size, n_cells,
                                                call = NULL) {
  msg <- sprintf(
    "All %d cells below minimum size %d (max observed: %d). Output panel would be empty.",
    n_cells, min_cell_size, max_observed_size
  )
  lwdid_error(msg,
              class = c("lwdid_insufficient_cell_size",
                        "lwdid_aggregation_error"),
              min_cell_size = min_cell_size,
              max_observed_size = max_observed_size,
              n_cells = n_cells, call = call)
}

# ============================================================================
# Warning Condition Convenience Constructors (5 types)
# ============================================================================

#' @rdname lwdid-conditions
#' @export
lwdid_small_sample_warning <- function(n, n_treated, n_control,
                                        cohort = NULL, period = NULL,
                                        call = NULL) {
  msg <- sprintf("Small sample size (N=%d, treated=%d, control=%d). Consider checking normality assumptions.", n, n_treated, n_control)
  if (!is.null(cohort)) {
    msg <- paste0(msg, sprintf(" Cohort=%s, period=%s.", as.character(cohort), as.character(period)))
  }
  lwdid_warning(msg, class = "lwdid_small_sample",
                n = n, n_treated = n_treated, n_control = n_control,
                cohort = cohort, period = period, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_overlap_warning <- function(n_extreme, ps_range, max_weight,
                                   call = NULL) {
  msg <- sprintf(
    "Overlap concern: %d units with extreme propensity scores (range: [%.4f, %.4f], max weight: %.2f).",
    n_extreme, ps_range[1], ps_range[2], max_weight
  )
  lwdid_warning(msg, class = "lwdid_overlap",
                n_extreme = n_extreme, ps_range = ps_range,
                max_weight = max_weight, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_numerical_warning <- function(detail, source_function,
                                     condition_number = NULL, call = NULL) {
  msg <- sprintf("Numerical issue in %s: %s", source_function, detail)
  if (!is.null(condition_number)) {
    msg <- paste0(msg, sprintf(" (condition number: %.2e)", condition_number))
  }
  lwdid_warning(msg, class = "lwdid_numerical",
                detail = detail, source_function = source_function,
                condition_number = condition_number, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_data_warning <- function(detail, action_taken, call = NULL) {
  msg <- sprintf("Data issue: %s. Action taken: %s", detail, action_taken)
  lwdid_warning(msg, class = "lwdid_data",
                detail = detail, action_taken = action_taken, call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_convergence_warning <- function(model, iterations, max_iterations,
                                       tolerance, call = NULL) {
  msg <- sprintf(
    "Convergence failure in %s: %d/%d iterations completed (tolerance: %.2e).",
    model, iterations, max_iterations, tolerance
  )
  lwdid_warning(msg, class = "lwdid_convergence",
                model = model, iterations = iterations,
                max_iterations = max_iterations, tolerance = tolerance,
                call = call)
}

#' @rdname lwdid-conditions
#' @export
lwdid_control_group_switch_warning <- function(original, switched_to,
                                                aggregate, call = NULL) {
  msg <- sprintf(
    paste0("Control group auto-switched from '%s' to 'never_treated' ",
           "for aggregate='%s'. Cohort/overall aggregation requires ",
           "never-treated units as the control group."),
    original, aggregate
  )
  lwdid_warning(msg, class = "lwdid_control_group_switch",
                original = original, switched_to = switched_to,
                aggregate = aggregate, call = call)
}
