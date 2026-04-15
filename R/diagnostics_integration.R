# ============================================================================
# diagnostics_integration.R -- Diagnostics extraction helpers and public API
#
# Story E8-07 starts by normalizing the diagnostics surface exposed by
# lwdid_result objects. Mainline producer wiring and standalone suites can
# layer on top of this file in later slices.
# ============================================================================

#' @title Extract diagnostics from lwdid results
#' @description Normalize and extract diagnostics attached to a
#'   \code{lwdid_result} object.
#' @param x A \code{lwdid_result} object.
#' @param type Character scalar indicating which diagnostics to return.
#'   One of \code{"all"}, \code{"clustering"}, \code{"selection"},
#'   \code{"trends"}, or \code{"sensitivity"}.
#' @return A normalized diagnostics list, a single diagnostics object, or
#'   \code{NULL} when diagnostics have not been run.
#' @export
get_diagnostics <- function(x, type = "all") {
  valid_types <- c("all", "clustering", "selection", "trends", "sensitivity")

  if (!(inherits(x, "lwdid_result") || inherits(x, "lwdid"))) {
    stop("x must inherit from 'lwdid_result' or 'lwdid'.", call. = FALSE)
  }

  if (!is.character(type) || length(type) != 1L || !type %in% valid_types) {
    stop(
      sprintf(
        "type must be one of: %s",
        paste(valid_types, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  diagnostics <- .normalize_integrated_diagnostics(x)
  has_diagnostics <- any(
    vapply(diagnostics, function(entry) !is.null(entry), logical(1))
  )

  if (!has_diagnostics) {
    message("Diagnostics were not run. Use lwdid(..., return_diagnostics=TRUE).")
    return(invisible(NULL))
  }

  if (identical(type, "all")) {
    return(diagnostics)
  }

  diagnostics[[type]]
}

#' @keywords internal
.attach_mainline_diagnostics <- function(result, validated) {
  params <- validated$validated_params
  never_treated_values <- validated$never_treated_values %||%
    params$never_treated_values %||%
    c(0, Inf)

  diagnostics <- result$diagnostics
  if (is.null(diagnostics)) {
    diagnostics <- list()
  }

  diagnostics$selection <- .run_suite_diagnostic("Selection diagnostics", function() {
    diagnose_selection_mechanism(
      validated$data,
      ivar = params$ivar,
      tvar = params$tvar,
      y = params$depvar,
      gvar = params$gvar,
      d = params$d,
      covariates = params$controls,
      never_treated_values = never_treated_values,
      verbose = FALSE
    )
  })

  diagnostics$trends <- .run_suite_diagnostic("Trend diagnostics", function() {
    lwdid_recommend_transformation(
      validated$data,
      y = params$depvar,
      ivar = params$ivar,
      tvar = params$tvar,
      gvar = params$gvar,
      controls = params$controls,
      never_treated_values = never_treated_values,
      run_all_diagnostics = TRUE,
      verbose = FALSE
    )
  })
  diagnostics$parallel_trends <- diagnostics$trends

  if (!is.null(params$gvar) && !is.null(params$cluster_var)) {
    diagnostics$clustering <- .run_suite_diagnostic("Clustering diagnostics", function() {
      diagnose_clustering(
        validated$data,
        ivar = params$ivar,
        potential_cluster_vars = c(params$cluster_var),
        gvar = params$gvar,
        d = params$d,
        verbose = FALSE
      )
    })
  }

  result$diagnostics <- diagnostics
  result
}

#' @keywords internal
.normalize_integrated_diagnostics <- function(x) {
  raw_diagnostics <- x$diagnostics
  if (is.null(raw_diagnostics)) {
    raw_diagnostics <- list()
  }

  trend_diag <- raw_diagnostics$trends
  if (is.null(trend_diag) && !is.null(raw_diagnostics$parallel_trends)) {
    trend_diag <- raw_diagnostics$parallel_trends
  }

  sensitivity_diag <- raw_diagnostics$sensitivity
  if (is.null(sensitivity_diag) && !is.null(x$sensitivity)) {
    sensitivity_diag <- x$sensitivity
  }

  list(
    clustering = raw_diagnostics$clustering %||% NULL,
    selection = raw_diagnostics$selection %||% NULL,
    trends = trend_diag %||% NULL,
    sensitivity = sensitivity_diag %||% NULL
  )
}

#' @title Run standalone diagnostics suite
#' @description Execute selection, trend, and optional clustering diagnostics
#'   without first fitting \code{lwdid()}.
#' @param data A panel \code{data.frame}.
#' @param y Outcome column name.
#' @param ivar Unit identifier column name.
#' @param tvar Time column name.
#' @param gvar Optional cohort column for staggered designs.
#' @param d Optional treatment indicator for common-timing designs.
#' @param controls Optional character vector of control-variable names.
#' @param cluster_vars Optional character vector of candidate cluster variables.
#'   When \code{NULL}, clustering diagnostics are skipped.
#' @param never_treated_values Numeric markers for never-treated cohorts.
#' @param verbose Logical scalar forwarded to the underlying diagnostics.
#' @return A \code{lwdid_diagnostics_suite} object.
#' @export
lwdid_diagnose <- function(
    data,
    y,
    ivar,
    tvar,
    gvar = NULL,
    d = NULL,
    controls = NULL,
    cluster_vars = NULL,
    never_treated_values = c(0, Inf),
    verbose = TRUE
) {
  # Always suppress sub-component printing during suite execution;

  results <- list(
    selection = .run_suite_diagnostic("Selection diagnostics", function() {
      diagnose_selection_mechanism(
        data,
        ivar = ivar,
        tvar = tvar,
        y = y,
        gvar = gvar,
        d = d,
        covariates = controls,
        never_treated_values = never_treated_values,
        verbose = FALSE
      )
    }),
    trends = .run_suite_diagnostic("Trend diagnostics", function() {
      lwdid_recommend_transformation(
        data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        controls = controls,
        never_treated_values = never_treated_values,
        run_all_diagnostics = TRUE,
        verbose = FALSE
      )
    })
  )

  if (!is.null(cluster_vars)) {
    results$clustering <- .run_suite_diagnostic("Clustering diagnostics", function() {
      diagnose_clustering(
        data,
        ivar = ivar,
        potential_cluster_vars = cluster_vars,
        gvar = gvar,
        d = d,
        verbose = FALSE
      )
    })
  }

  structure(results, class = c("lwdid_diagnostics_suite", "list"))
}

#' Print method for lwdid diagnostics suite
#'
#' @param x An object of class \code{lwdid_diagnostics_suite}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.lwdid_diagnostics_suite <- function(x, ...) {
  cat("lwdid diagnostics suite\n")
  cat(strrep("=", 50), "\n")

  selection_risk <- if (!is.null(x$selection)) {
    x$selection$selection_risk %||% "unknown"
  } else {
    "Not run"
  }
  cat(sprintf("Selection diagnostics: %s\n", selection_risk))

  trend_summary <- if (!is.null(x$trends)) {
    sprintf(
      "%s (%s)",
      x$trends$recommended_method %||% "unknown",
      x$trends$confidence_level %||% "unknown"
    )
  } else {
    "Not run"
  }
  cat(sprintf("Trend diagnostics: %s\n", trend_summary))

  clustering_summary <- if (!is.null(x$clustering)) {
    x$clustering$recommended_cluster_var %||% "none"
  } else {
    "Not run"
  }
  cat(sprintf("Clustering diagnostics: %s\n", clustering_summary))

  invisible(x)
}

#' Summary method for lwdid diagnostics suite
#'
#' @param object An object of class \code{lwdid_diagnostics_suite}.
#' @param ... Additional arguments passed to component summaries.
#' @return Invisibly returns \code{object}.
#' @export
summary.lwdid_diagnostics_suite <- function(object, ...) {
  cat("Diagnostics suite summary\n")
  cat(strrep("=", 50), "\n\n")

  .emit_suite_component_summary("Selection diagnostics", object$selection, ...)
  .emit_suite_component_summary("Trend diagnostics", object$trends, ...)
  .emit_suite_component_summary("Clustering diagnostics", object$clustering, ...)

  invisible(object)
}

#' @keywords internal
.run_suite_diagnostic <- function(label, diagnostic_call) {
  tryCatch(
    diagnostic_call(),
    error = function(e) {
      warning(
        sprintf("%s failed: %s", label, conditionMessage(e)),
        call. = FALSE
      )
      NULL
    }
  )
}

#' @keywords internal
.emit_suite_component_summary <- function(label, component, ...) {
  cat(label, "\n", sep = "")
  cat(strrep("-", 20), "\n", sep = "")

  if (is.null(component)) {
    cat("Not run\n\n")
    return(invisible(NULL))
  }

  rendered <- utils::capture.output(summary(component, ...))
  if (length(rendered) == 0L) {
    rendered <- utils::capture.output(print(component, ...))
  }
  cat(paste(rendered, collapse = "\n"), "\n\n", sep = "")

  invisible(NULL)
}
