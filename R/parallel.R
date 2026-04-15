# parallel.R
# Performance optimization and parallel computation support for lwdid.
# Provides future-based parallel execution for staggered DiD estimation,
# RI permutation inference, and Wild Cluster Bootstrap.
#
# Architecture:
#   setup_parallel_internal() - configure future plan (internal)
#   run_parallel()            - unified parallel loop executor
#   safe_iterate()            - error/warning-capturing wrapper
#   setup_parallel()          - user-facing plan configuration
#   .can_use_parallel()       - check if parallel backend active


# =========================================================================
# safe_iterate: error/warning-capturing single-iteration wrapper
# =========================================================================

#' Safely execute a single iteration with error and warning capture
#'
#' Wraps a function call to capture errors and warnings, returning a
#' structured result. Used by \code{\link{run_parallel}} for fault
#' tolerance in parallel loops.
#'
#' @param fun Function to execute (zero-argument closure).
#' @return A list with elements:
#'   \describe{
#'     \item{result}{The return value of \code{fun} (NULL on failure).}
#'     \item{status}{Character: \code{"success"} or \code{"failed"}.}
#'     \item{error}{Character error message (NULL on success).}
#'     \item{warnings}{Character vector of captured warning messages.}
#'   }
#' @keywords internal
safe_iterate <- function(fun) {
  warnings_collected <- character(0L)

  result <- tryCatch(
    withCallingHandlers(
      list(result = fun(), status = "success",
           error = NULL),
      warning = function(w) {
        warnings_collected <<- c(warnings_collected,
                                 conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(result = NULL, status = "failed",
           error = conditionMessage(e))
    }
  )

  result$warnings <- warnings_collected
  result
}


# =========================================================================
# setup_parallel_internal: internal plan configuration
# =========================================================================

#' Configure parallel backend for internal use
#'
#' Sets up the \pkg{future} plan based on the \code{parallel} and
#' \code{n_cores} arguments. Returns a cleanup function that restores
#' the original plan, suitable for use with \code{on.exit()}.
#'
#' @param parallel Logical. Whether to enable parallel execution.
#' @param n_cores Integer or NULL. Number of workers. NULL means auto-detect.
#' @return A function that restores the original plan (NULL if no change).
#' @keywords internal
.setup_parallel_internal <- function(parallel = FALSE,
                                     n_cores = NULL) {
  if (!isTRUE(parallel)) {
    return(invisible(NULL))
  }

  # Check required packages
  pkgs <- c("future", "future.apply")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace,
                               logical(1L), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    stop_lwdid(
      sprintf(
        "Parallel computation requires: %s. Install with: install.packages(c(%s))",
        paste(missing_pkgs, collapse = ", "),
        paste0('"', missing_pkgs, '"', collapse = ", ")
      ),
      class = "lwdid_input"
    )
  }

  # If a non-sequential plan is already active, respect it
  current_plan <- tryCatch(future::plan(), error = function(...) NULL)
  if (!is.null(current_plan) &&
      !inherits(current_plan, "sequential")) {
    return(invisible(NULL))
  }

  # Auto-detect cores
  if (is.null(n_cores)) {
    n_cores <- max(1L, future::availableCores() - 1L)
  }
  n_cores <- as.integer(n_cores)
  if (n_cores < 1L) {
    stop_lwdid(
      sprintf("n_cores must be >= 1, got: %d", n_cores),
      class = "lwdid_input"
    )
  }

  # CRAN check environment: limit to 2 cores
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    n_cores <- min(n_cores, 2L)
  }

  if (n_cores == 1L) {
    return(invisible(NULL))
  }

  # Save current plan stack for restoration
  old_plan <- future::plan("list")

  # Auto-select strategy
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, workers = n_cores)
  } else if (interactive() && nzchar(Sys.getenv("RSTUDIO"))) {
    future::plan(future::multisession, workers = n_cores)
  } else {
    future::plan(future::multicore, workers = n_cores)
  }

  # Return cleanup function
  function() future::plan(old_plan)
}


# =========================================================================
# run_parallel: unified parallel loop executor
# =========================================================================

#' Execute a function over a sequence with optional parallelism
#'
#' Unified parallel/sequential loop executor with fault tolerance.
#' When \code{parallel = TRUE} and \pkg{future.apply} is available,
#' uses \code{future.apply::future_lapply}. Otherwise, falls back to
#' \code{base::lapply}. Individual iteration failures are captured and
#' reported without terminating the entire loop.
#'
#' @param X Vector or list to iterate over.
#' @param FUN Function to apply to each element of \code{X}.
#' @param ... Additional arguments passed to every call of \code{FUN}.
#' @param parallel Logical. Whether to use parallel execution.
#' @param n_cores Integer or NULL. Number of parallel workers.
#' @param fail_threshold Numeric. Maximum fraction of failed iterations
#'   before issuing a warning (default 0.1).
#' @param future.seed Logical or integer. Seed control for
#'   \code{future_lapply} (default TRUE for reproducibility).
#' @param future.scheduling Numeric. Task chunking strategy
#'   (default 1.0; use 2.0 for uneven workloads).
#' @param task_name Character. Label for progress/error messages.
#'
#' @return A list of results from successful iterations.
#'
#' @details
#' Iterations that throw errors are captured by \code{\link{safe_iterate}}
#' and excluded from the result. If more than 50\% of iterations fail,
#' an error is raised. If more than \code{fail_threshold} fail, a warning
#' is issued.
#'
#' @examples
#' # Sequential execution
#' run_parallel(1:5, function(x) x^2, parallel = FALSE)
#'
#' @seealso \code{\link{setup_parallel}}, \code{\link{safe_iterate}}
#' @export
run_parallel <- function(X, FUN, ...,
                         parallel = FALSE,
                         n_cores = NULL,
                         fail_threshold = 0.1,
                         future.seed = TRUE,
                         future.scheduling = 1.0,
                         task_name = "iteration") {
  n_total <- length(X)

  # Empty input fast return
  if (n_total == 0L) return(list())

  # Set up parallel environment
  cleanup <- .setup_parallel_internal(
    parallel = parallel, n_cores = n_cores
  )
  if (!is.null(cleanup)) on.exit(cleanup(), add = TRUE)

  # Capture extra args to avoid ... propagation issues in closures
  extra_args <- list(...)

  # Wrap FUN with safe_iterate for fault tolerance
  safe_FUN <- function(x) {
    safe_iterate(function() do.call(FUN, c(list(x), extra_args)))
  }

  # Determine actual execution mode
  use_par <- isTRUE(parallel) && .can_use_parallel()

  if (use_par) {
    results <- future.apply::future_lapply(
      X, safe_FUN,
      future.seed = future.seed,
      future.scheduling = future.scheduling
    )
  } else {
    results <- lapply(X, safe_FUN)
  }

  # Tally failures
  statuses <- vapply(results, `[[`, character(1L), "status")
  n_failed <- sum(statuses == "failed")

  if (n_failed > 0L) {
    fail_ratio <- n_failed / n_total
    error_msgs <- vapply(
      results[statuses == "failed"],
      function(r) r$error, character(1L)
    )
    unique_errors <- unique(error_msgs)

    if (fail_ratio > 0.5) {
      stop_lwdid(
        sprintf(
          "%s: %d/%d iterations failed (%.0f%%), exceeds 50%% threshold",
          task_name, n_failed, n_total, fail_ratio * 100
        ),
        class = "lwdid_input"
      )
    }

    if (fail_ratio > fail_threshold) {
      warn_lwdid(
        sprintf(
          "%s: %d/%d iterations failed (%.0f%%)",
          task_name, n_failed, n_total, fail_ratio * 100
        ),
        class = "lwdid_data"
      )
    }
  }

  # Extract successful results
  successful <- results[statuses == "success"]
  lapply(successful, `[[`, "result")
}


# =========================================================================
# User-facing setup_parallel
# =========================================================================

#' Set up parallel computation plan
#'
#' Configures the parallel backend using the \pkg{future} package.
#' Call this function before running \code{\link{lwdid}} with
#' \code{parallel = TRUE} to control the number of workers and strategy.
#'
#' @param n_workers Integer or NULL. Number of parallel workers. If NULL,
#'   defaults to \code{future::availableCores() - 1} (at least 1).
#' @param strategy Character(1). Parallelization strategy passed to
#'   \code{future::plan()}. One of \code{"multisession"} (default,
#'   works on all platforms), \code{"multicore"} (Unix/macOS only,
#'   faster but not safe in RStudio/GUI), or \code{"sequential"}
#'   (disable parallelism).
#'
#' @return Invisibly returns the previous \code{future::plan()}.
#'
#' @details
#' This function requires the \pkg{future} package to be installed.
#' If not installed, an informative error is raised.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("future", quietly = TRUE)) {
#'   setup_parallel(n_workers = 2)
#'   # Now lwdid(..., parallel = TRUE) will use 2 workers
#'   # Reset to sequential when done:
#'   setup_parallel(strategy = "sequential")
#' }
#' }
#'
#' @seealso \code{\link{lwdid}}, \code{\link{run_parallel}}
#' @export
setup_parallel <- function(n_workers = NULL,
                           strategy = c("multisession", "multicore",
                                        "sequential")) {
  if (!requireNamespace("future", quietly = TRUE)) {
    stop(
      "Package 'future' is required for parallel computation.\n",
      "Install it with: install.packages('future')",
      call. = FALSE
    )
  }

  strategy <- match.arg(strategy)

  if (identical(strategy, "sequential")) {
    old_plan <- future::plan(future::sequential)
    message("lwdid: parallel computation disabled (sequential mode).")
    return(invisible(old_plan))
  }

  if (is.null(n_workers)) {
    n_workers <- max(1L, future::availableCores() - 1L)
  }
  n_workers <- as.integer(n_workers)
  if (n_workers < 1L) n_workers <- 1L

  # CRAN check environment: limit to 2 cores
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    n_workers <- min(n_workers, 2L)
  }

  plan_fn <- switch(
    strategy,
    multisession = future::multisession,
    multicore    = future::multicore,
    future::multisession
  )

  old_plan <- future::plan(plan_fn, workers = n_workers)
  message(
    sprintf(
      "lwdid: parallel computation enabled (%s, %d workers).",
      strategy, n_workers
    )
  )
  invisible(old_plan)
}


# =========================================================================
# Internal helpers
# =========================================================================

#' Check whether parallel execution is available and active
#' @return Logical scalar.
#' @keywords internal
.can_use_parallel <- function() {
  if (!requireNamespace("future", quietly = TRUE)) return(FALSE)
  if (!requireNamespace("future.apply", quietly = TRUE)) return(FALSE)

  current_plan <- tryCatch(
    future::plan(),
    error = function(...) NULL
  )
  if (is.null(current_plan)) return(FALSE)

  !inherits(current_plan, "sequential")
}


#' Parallel-aware lapply with graceful fallback
#'
#' Internal utility that applies a function over a list, using parallel
#' execution when available and falling back to sequential lapply otherwise.
#'
#' @param X List or vector to iterate over.
#' @param FUN Function to apply.
#' @param ... Additional arguments to \code{FUN}.
#' @return A list of results.
#' @keywords internal
.parallel_lapply <- function(X, FUN, ...) {
  if (.can_use_parallel()) {
    future.apply::future_lapply(X, FUN, ..., future.seed = TRUE)
  } else {
    lapply(X, FUN, ...)
  }
}
