#' Resolve the effective seasonal variable name
#'
#' @param season_var Optional seasonal variable.
#' @param quarter Optional backward-compatible alias.
#' @return Character scalar naming the seasonal column.
#' @keywords internal
.resolve_season_var <- function(season_var = NULL, quarter = NULL) {
  effective <- season_var
  if (is.null(effective)) {
    effective <- quarter
  }

  if (is.null(effective) || !is.character(effective) || length(effective) != 1L) {
    stop_lwdid(
      message = paste(
        "Seasonal transformations require 'season_var' or the backward-compatible",
        "'quarter' parameter."
      ),
      class = "lwdid_invalid_parameter",
      param = "season_var",
      value = effective,
      allowed = "character(1)"
    )
  }

  effective
}

#' Apply unit-level pre-period exclusion for seasonal validation
#'
#' @param data data.frame-like object.
#' @param ivar Unit identifier column name.
#' @param tindex Time index column name.
#' @param post Post indicator column name.
#' @param exclude_pre_periods Integer number of trailing pre-periods to exclude.
#' @return Integer vector aligned with `data[[post]]`.
#' @keywords internal
.compute_effective_seasonal_post <- function(
  data, ivar, tindex, post, exclude_pre_periods = 0L
) {
  effective_post <- as.integer(data[[post]])

  if (exclude_pre_periods <= 0L) {
    return(effective_post)
  }

  for (uid in unique(data[[ivar]])) {
    unit_mask <- data[[ivar]] == uid
    unit_pre_mask <- unit_mask & effective_post == 0L & !is.na(data[[tindex]])
    pre_times <- sort(unique(data[[tindex]][unit_pre_mask]))

    if (length(pre_times) > exclude_pre_periods) {
      times_to_exclude <- utils::tail(pre_times, exclude_pre_periods)
      exclude_mask <- unit_mask & data[[tindex]] %in% times_to_exclude
      effective_post[exclude_mask] <- 1L
    }
  }

  effective_post
}

#' Validate seasonal values are encoded as integers in 1:Q
#'
#' @param values Seasonal values.
#' @param season_var Seasonal column name.
#' @param Q Number of seasonal buckets.
#' @return Integer vector of validated seasonal values.
#' @keywords internal
.normalize_season_values <- function(values, season_var, Q) {
  raw_values <- unique(values[!is.na(values)])
  if (length(raw_values) == 0L) {
    return(integer(0))
  }

  numeric_values <- suppressWarnings(as.numeric(raw_values))
  invalid_numeric <- is.na(numeric_values) |
    abs(numeric_values - round(numeric_values)) > sqrt(.Machine$double.eps)
  integer_values <- as.integer(round(numeric_values))
  invalid_range <- !invalid_numeric & (integer_values < 1L | integer_values > Q)

  if (any(invalid_numeric) || any(invalid_range)) {
    invalid_values <- raw_values[invalid_numeric | invalid_range]
    hint <- if (any(!invalid_numeric & integer_values == 0L)) {
      paste0(
        " If your seasonal values are encoded as 0-", Q - 1L,
        ", recode them to 1-", Q, " before calling lwdid()."
      )
    } else {
      ""
    }

    stop_lwdid(
      message = paste0(
        "Seasonal column '", season_var, "' contains invalid values: ",
        paste(sort(unique(as.character(invalid_values))), collapse = ", "),
        ". Expected integer values in {1, 2, ..., ", Q, "}.",
        hint
      ),
      class = "lwdid_invalid_parameter",
      param = "season_var",
      value = invalid_values,
      allowed = seq_len(Q)
    )
  }

  integer_values
}

#' Compute seasonal parameter requirements for validation and transforms
#'
#' @param method Seasonal rolling method.
#' @param n_unique_seasons Integer count of observed pre-treatment seasons.
#' @return Named list with parameter-count contract.
#' @keywords internal
.seasonal_parameter_requirements <- function(method, n_unique_seasons) {
  if (!method %in% c("demeanq", "detrendq")) {
    stop_lwdid(
      message = sprintf(
        "Seasonal parameter counting only supports demeanq/detrendq. Got '%s'.",
        as.character(method)
      ),
      class = "lwdid_invalid_parameter",
      param = "method",
      value = method,
      allowed = c("demeanq", "detrendq")
    )
  }

  n_unique_seasons <- as.integer(n_unique_seasons)
  if (length(n_unique_seasons) != 1L || is.na(n_unique_seasons) || n_unique_seasons < 0L) {
    stop_lwdid(
      message = sprintf(
        "n_unique_seasons must be a non-negative integer. Got: %s",
        paste(as.character(n_unique_seasons), collapse = ", ")
      ),
      class = "lwdid_invalid_parameter",
      param = "n_unique_seasons",
      value = n_unique_seasons,
      allowed = "non-negative integer(1)"
    )
  }

  n_params <- n_unique_seasons + if (identical(method, "detrendq")) 1L else 0L

  list(
    method = method,
    n_unique_seasons = n_unique_seasons,
    n_params = as.integer(n_params),
    min_required = as.integer(n_params + 1L)
  )
}

#' Check whether post-treatment seasons are covered in the pre-period
#'
#' @param data data.frame-like object.
#' @param season_var Seasonal column name.
#' @param ivar Unit identifier column name.
#' @param post Post indicator column name.
#' @param Q Number of seasonal buckets.
#' @return Structured coverage summary.
#' @keywords internal
.check_season_coverage <- function(data, season_var, ivar, post, Q = 4L,
                                   registry = NULL) {
  freq_label <- FREQ_LABELS[[as.character(Q)]]
  if (is.null(freq_label)) {
    freq_label <- "season"
  }

  units_with_gaps <- list()
  total_uncovered <- 0L
  warnings_seen <- character(0)

  for (uid in unique(data[[ivar]])) {
    unit_mask <- data[[ivar]] == uid
    pre_seasons <- unique(data[[season_var]][unit_mask & data[[post]] == 0L])
    pre_seasons <- sort(pre_seasons[!is.na(pre_seasons)])
    post_seasons <- unique(data[[season_var]][unit_mask & data[[post]] == 1L])
    post_seasons <- sort(post_seasons[!is.na(post_seasons)])
    uncovered <- setdiff(post_seasons, pre_seasons)

    if (length(uncovered) > 0L) {
      units_with_gaps[[as.character(uid)]] <- unname(uncovered)
      total_uncovered <- total_uncovered + length(uncovered)

      warning_message <- sprintf(
        paste(
          "Unit %s has post-treatment %s(s) %s that were not observed in the",
          "pre-treatment window %s. R records this as a warning and allows",
          "the seasonal gate to pass."
        ),
        as.character(uid),
        freq_label,
        paste(uncovered, collapse = ", "),
        paste(pre_seasons, collapse = ", ")
      )

      warnings_seen <- c(warnings_seen, warning_message)
      if (!is.null(registry) && is.function(registry$register)) {
        registry$register(
          category = "lwdid_data",
          message = warning_message,
          context = list(
            detail = "seasonal_coverage_gap",
            unit = as.character(uid),
            uncovered = paste(uncovered, collapse = ", "),
            pre_window = paste(pre_seasons, collapse = ", ")
          )
        )
      } else {
        warn_lwdid(
          message = warning_message,
          class = "lwdid_data",
          detail = "seasonal_coverage_gap",
          action_taken = "coverage recorded as warning"
        )
      }
    }
  }

  list(
    all_covered = length(units_with_gaps) == 0L,
    units_with_gaps = units_with_gaps,
    total_uncovered = as.integer(total_uncovered),
    warnings = unname(unique(warnings_seen))
  )
}

#' Auto-detect seasonal frequency for story-level seasonal validation
#'
#' @param data data.frame-like object.
#' @param tvar Time variable column name.
#' @param ivar Optional unit identifier column name.
#' @return Named list with `Q`, `frequency`, and `confidence`, or NULL.
#' @keywords internal
.auto_detect_frequency <- function(data, tvar, ivar = NULL) {
  dt <- if (data.table::is.data.table(data)) {
    data
  } else {
    data.table::as.data.table(data)
  }

  detection <- .detect_frequency(dt, tvar, ivar)
  if (is.null(detection$frequency) || is.null(detection$Q)) {
    return(NULL)
  }

  list(
    Q = as.integer(detection$Q),
    frequency = detection$frequency,
    confidence = as.numeric(detection$confidence)
  )
}

#' Validate inputs for seasonal rolling transformations
#'
#' @param data data.frame-like object.
#' @param y Outcome column name.
#' @param season_var Optional seasonal column name.
#' @param Q Number of seasonal buckets.
#' @param ivar Unit identifier column name.
#' @param tindex Time index column name.
#' @param post Post indicator column name.
#' @param method Seasonal rolling method.
#' @param exclude_pre_periods Integer number of pre-periods to exclude.
#' @param min_global_pre_periods Minimum unique pre-treatment periods required.
#' @param quarter Optional backward-compatible alias for season_var.
#' @return Structured seasonal validation summary.
#' @keywords internal
.validate_seasonal_inputs <- function(
  data, y, season_var = NULL, Q, ivar, tindex, post = NULL,
  method = "demeanq", exclude_pre_periods = 0L,
  min_global_pre_periods = NULL, quarter = NULL, registry = NULL
) {
  .validate_positive_integer(Q, "Q", min_val = 2L)
  effective_season_var <- .resolve_season_var(season_var = season_var, quarter = quarter)
  min_global_pre_periods <- min_global_pre_periods %||%
    if (identical(method, "detrendq")) 2L else 1L

  if (!effective_season_var %in% names(data)) {
    stop_lwdid(
      message = sprintf(
        "Required column '%s' not found in data. Available columns: %s",
        effective_season_var,
        paste(names(data), collapse = ", ")
      ),
      class = "lwdid_missing_column",
      column = effective_season_var,
      available = names(data)
    )
  }

  if (is.null(post) || !post %in% names(data)) {
    stop_lwdid(
      message = "Seasonal validation requires a valid post indicator column.",
      class = "lwdid_invalid_parameter",
      param = "post",
      value = post,
      allowed = "existing column name"
    )
  }

  .normalize_season_values(data[[effective_season_var]], effective_season_var, Q)

  validated_data <- as.data.frame(data, stringsAsFactors = FALSE)
  validated_data[[".seasonal_post"]] <- .compute_effective_seasonal_post(
    data = validated_data,
    ivar = ivar,
    tindex = tindex,
    post = post,
    exclude_pre_periods = as.integer(exclude_pre_periods)
  )

  pre_rows <- validated_data[
    validated_data[[".seasonal_post"]] == 0L,
    ,
    drop = FALSE
  ]

  if (nrow(pre_rows) > 0L && all(is.na(pre_rows[[tindex]]))) {
    stop_lwdid(
      message = sprintf(
        paste(
          "All pre-treatment tindex values are missing.",
          "rolling('%s') requires valid time index values."
        ),
        method
      ),
      class = c("lwdid_insufficient_pre_periods", "lwdid_insufficient_data"),
      rolling = method
    )
  }

  pre_tindex <- pre_rows[[tindex]][!is.na(pre_rows[[tindex]])]
  n_pre_periods <- length(unique(pre_tindex))

  if (n_pre_periods < min_global_pre_periods) {
    stop_lwdid(
      message = sprintf(
        "rolling('%s') requires at least %d pre-treatment period(s). Found: %d.",
        method,
        as.integer(min_global_pre_periods),
        as.integer(n_pre_periods)
      ),
      class = c("lwdid_insufficient_pre_periods", "lwdid_insufficient_data"),
      available = as.integer(n_pre_periods),
      required = as.integer(min_global_pre_periods),
      rolling = method
    )
  }

  param_offset <- if (identical(method, "demeanq")) 0L else 1L
  valid_units <- 0L

  for (uid in unique(validated_data[[ivar]])) {
    unit_mask <- validated_data[[ivar]] == uid
    unit_pre <- validated_data[unit_mask & validated_data[[".seasonal_post"]] == 0L, , drop = FALSE]
    unit_pre_count <- nrow(unit_pre)

    if (unit_pre_count < min_global_pre_periods) {
      stop_lwdid(
        message = sprintf(
          "Unit %s has only %d pre-period observation(s). rolling('%s') requires at least %d.",
          as.character(uid),
          unit_pre_count,
          method,
          as.integer(min_global_pre_periods)
        ),
        class = c("lwdid_insufficient_pre_periods", "lwdid_insufficient_data"),
        unit = uid,
        available = unit_pre_count,
        required = as.integer(min_global_pre_periods),
        rolling = method
      )
    }

    valid_mask <- !is.na(unit_pre[[y]]) & !is.na(unit_pre[[effective_season_var]])
    if (identical(method, "detrendq")) {
      valid_mask <- valid_mask & !is.na(unit_pre[[tindex]])
    }

    requirements <- .seasonal_parameter_requirements(
      method = method,
      n_unique_seasons = length(unique(unit_pre[[effective_season_var]][valid_mask]))
    )
    n_unique_seasons <- requirements$n_unique_seasons
    min_required <- requirements$min_required

    if (unit_pre_count < min_required) {
      stop_lwdid(
        message = sprintf(
          paste(
            "Unit %s has %d pre-period observation(s) with %d distinct season(s).",
            "rolling('%s') requires at least %d observation(s) so that df >= 1."
          ),
          as.character(uid),
          unit_pre_count,
          n_unique_seasons,
          method,
          min_required
        ),
        class = c("lwdid_insufficient_pre_periods", "lwdid_insufficient_data"),
        unit = uid,
        available = unit_pre_count,
        required = min_required,
        rolling = method
      )
    }

    valid_units <- valid_units + 1L
  }

  coverage <- .check_season_coverage(
    data = validated_data,
    season_var = effective_season_var,
    ivar = ivar,
    post = ".seasonal_post",
    Q = Q,
    registry = registry
  )

  if (length(pre_tindex) == 0L) {
    stop_lwdid(
      message = "Cannot compute K because all pre-treatment tindex values are missing.",
      class = c("lwdid_insufficient_pre_periods", "lwdid_insufficient_data"),
      rolling = method
    )
  }

  list(
    is_valid = TRUE,
    K = as.integer(max(pre_tindex)),
    n_units_valid = as.integer(valid_units),
    n_units_invalid = 0L,
    invalid_units = character(0),
    season_coverage = coverage[c("all_covered", "units_with_gaps", "total_uncovered")],
    warnings = coverage$warnings
  )
}
