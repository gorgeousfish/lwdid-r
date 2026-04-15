# ===========================================================================
# zzz.R — Package load hooks
#
# Implements .onLoad() and .onAttach() hooks.
# .onLoad(): Initializes package-level mutable state (global WarningRegistry,
#            version cache, package options)
# .onAttach(): Displays startup message
# ===========================================================================

#' @title Package Load Hooks
#' @description Package-level environment and lifecycle management for lwdid.
#' @name lwdid-zzz
#' @importFrom stats as.formula binomial coef complete.cases glm hatvalues lm
#'   lm.fit model.frame pnorm predict pt qnorm qt residuals rnorm runif sd
#'   setNames var vcov
#' @importFrom utils head tail
#' @keywords internal
NULL

# Suppress R CMD check NOTEs for data.table and ggplot2 NSE variables
utils::globalVariables(c(
  ".", ".N", ".SD", ":=",
  "..required_vars",
  "d_", "post_", "tindex", "tq",
  "is_nt", "n", "n_unique", "post_val", "gval",
  "g",
  "D_ig", "N", "alpha_hat", "beta_hat", "pre_mean",
  "ref_val", "y_trans", "y_trans_pre",
  ".data", "cohort",
  # plot_diagnostics.R NSE variables
  "n_pre", "att", "ci_lower", "ci_upper", "significant",
  "excluded", "fitted", "is_recommended", "dimension",
  "sub_score", "score", "method", "slope_label",
  "time", "is_control",
  # plot.R NSE variables (Story 10.5)
  "y", "group", "xint", "unit_id", "size",
  "cluster", "period", "coefficient", "mean_y",
  "attrition_rate", "label",
  # R CMD check --as-cran: remaining NSE variables
  ".trend_treat_dummy", "_gvar_dummy", "_missing", "_y_lag",
  "event_time", "i.ydot_postavg", "n_missing", "n_y",
  "post_placebo", "pvalue", "sd_y", "se_y",
  "t_crit", "treat_cohort", "treatment_time", "ydot_postavg",
  # preprocessing.R NSE variables
  "_n_obs", "_ess", "_outcome_var",
  "..group_cols", "..time_var", "sd_treat", "sd_gvar"
))

# Required for data.table to recognize this package as data.table-aware
.datatable.aware <- TRUE

# ---------------------------------------------------------------------------
# Package-level private environment (stores global mutable state)
# parent = emptyenv() ensures no inheritance from global environment
# Environment objects have reference semantics (no copy-on-modify)
# Slots (initialized by .onLoad()):
#   - warning_registry: global WarningRegistry singleton
#   - version: cached package version (package_version object)
# ---------------------------------------------------------------------------
.lwdid_env <- new.env(parent = emptyenv())

# ---------------------------------------------------------------------------
# .onLoad: Called when package namespace is loaded
# ---------------------------------------------------------------------------
.onLoad <- function(libname, pkgname) {
  # 1. Initialize global warning registry
  .lwdid_env$warning_registry <- new_warning_registry()

  # 2. Cache package version
  .lwdid_env$version <- tryCatch(
    utils::packageVersion(pkgname),
    error = function(e) package_version("0.0.0.9000")
  )

  # 3. Set default package options (guard mode, do not override user presets)
  #    Valid values: "quiet", "default", "verbose"
  if (is.null(getOption("lwdid.verbose"))) {
    options(lwdid.verbose = "default")
  }

  trend_s3_methods <- list(
    list("print", "lwdid_parallel_trends", print.lwdid_parallel_trends),
    list("summary", "lwdid_parallel_trends", summary.lwdid_parallel_trends),
    list("plot", "lwdid_parallel_trends", plot.lwdid_parallel_trends),
    list("print", "lwdid_heterogeneous_trends", print.lwdid_heterogeneous_trends),
    list("summary", "lwdid_heterogeneous_trends", summary.lwdid_heterogeneous_trends),
    list("plot", "lwdid_heterogeneous_trends", plot.lwdid_heterogeneous_trends),
    list("print", "lwdid_transformation_recommendation", print.lwdid_transformation_recommendation),
    list("summary", "lwdid_transformation_recommendation", summary.lwdid_transformation_recommendation),
    list("plot", "lwdid_transformation_recommendation", plot.lwdid_transformation_recommendation),
    list("print", "lwdid_diagnostics_suite", print.lwdid_diagnostics_suite),
    list("summary", "lwdid_diagnostics_suite", summary.lwdid_diagnostics_suite)
  )

  namespace_env <- asNamespace(pkgname)
  for (method_info in trend_s3_methods) {
    base::registerS3method(
      method_info[[1L]],
      method_info[[2L]],
      method_info[[3L]],
      envir = namespace_env
    )
  }

  invisible()
}

# ---------------------------------------------------------------------------
# .onAttach: Called when package is attached to search path (via library())
# ---------------------------------------------------------------------------
.onAttach <- function(libname, pkgname) {
  ver <- .lwdid_env$version  # use cached version
  packageStartupMessage(
    sprintf("lwdid %s: Lee-Wooldridge DiD Estimation", ver)
  )
}

# ---------------------------------------------------------------------------
# Internal helper: access the package-level private environment
# ---------------------------------------------------------------------------
#' Get the lwdid package-level private environment
#'
#' Returns the package-level private environment used for storing
#' the global WarningRegistry and cached version information.
#'
#' @return Environment with \code{warning_registry} and \code{version} slots.
#' @keywords internal
get_lwdid_env <- function() {
  .lwdid_env
}
