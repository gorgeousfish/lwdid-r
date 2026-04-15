# ============================================================================
# clustering_diagnostics.R — Internal clustering diagnostics primitives
#
# Foundations for clustering-structure diagnostics in Lee-Wooldridge DiD
# applications. This file currently provides Task E8-04.1 internals:
# input validation, cluster-level classification, and single-variable
# statistical analysis.
# ============================================================================

#' @title Clustering Diagnostics Internals
#' @description Internal helpers for clustering diagnostics.
#' @name clustering_diagnostics
#' @keywords internal
NULL


#' Validate inputs for clustering diagnostics
#'
#' @param data data.frame or data.table.
#' @param ivar character(1). Unit identifier.
#' @param potential_cluster_vars character vector of candidate cluster variables.
#' @param gvar character(1) or NULL. Cohort variable for staggered designs.
#' @param d character(1) or NULL. Treatment indicator for common timing.
#' @return NULL (invisible). Raises errors on failure.
#' @keywords internal
.validate_clustering_inputs <- function(
    data, ivar, potential_cluster_vars, gvar, d
) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }

  if (!is.character(ivar) || length(ivar) != 1L || !ivar %in% names(data)) {
    stop(sprintf("Unit variable '%s' not found in data", ivar), call. = FALSE)
  }

  if (!is.character(potential_cluster_vars) ||
      length(potential_cluster_vars) == 0L) {
    stop(
      "At least one potential cluster variable must be specified",
      call. = FALSE
    )
  }

  missing_clusters <- setdiff(potential_cluster_vars, names(data))
  if (length(missing_clusters) > 0L) {
    stop(
      sprintf(
        "Cluster variable '%s' not found in data",
        missing_clusters[[1L]]
      ),
      call. = FALSE
    )
  }

  if (!is.null(gvar) && (!is.character(gvar) || length(gvar) != 1L ||
                         !gvar %in% names(data))) {
    stop(sprintf("Cohort variable '%s' not found in data", gvar), call. = FALSE)
  }

  if (!is.null(d) && (!is.character(d) || length(d) != 1L ||
                      !d %in% names(data))) {
    stop(sprintf("Treatment variable '%s' not found in data", d), call. = FALSE)
  }

  if (is.null(gvar) && is.null(d)) {
    stop("Either gvar or d must be specified", call. = FALSE)
  }

  invisible(NULL)
}


#' Determine clustering level relative to unit
#'
#' @param data data.frame or data.table.
#' @param ivar character(1). Unit identifier.
#' @param cluster_var character(1). Candidate clustering variable.
#' @return character(1). One of "higher", "same", or "lower".
#' @keywords internal
.determine_cluster_level <- function(data, ivar, cluster_var) {
  dt <- data.table::as.data.table(data)
  n_units <- data.table::uniqueN(dt[[ivar]])
  n_clusters <- data.table::uniqueN(dt[[cluster_var]])
  units_per_cluster <- dt[, .(n_units = data.table::uniqueN(get(ivar))), by = cluster_var]

  if (n_clusters < n_units && any(units_per_cluster$n_units > 1L)) {
    return("higher")
  }
  if (n_clusters == n_units) {
    return("same")
  }
  "lower"
}


#' Analyze a single clustering variable
#'
#' @param data data.frame or data.table.
#' @param ivar character(1). Unit identifier.
#' @param cluster_var character(1). Candidate clustering variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param d character(1) or NULL. Treatment indicator.
#' @return Named list with Task E8-04.1 diagnostics.
#' @keywords internal
.analyze_cluster_var <- function(data, ivar, cluster_var, gvar, d) {
  dt <- data.table::as.data.table(data)

  cluster_sizes_dt <- dt[, .(N = .N), by = cluster_var]
  cluster_ids <- as.character(cluster_sizes_dt[[cluster_var]])
  cluster_sizes <- stats::setNames(as.integer(cluster_sizes_dt$N), cluster_ids)
  n_clusters <- as.integer(length(cluster_sizes))

  if (!is.null(gvar)) {
    treated_mask <- !is_never_treated(dt[[gvar]])
    treat_var <- gvar
  } else if (!is.null(d)) {
    treated_mask <- !is.na(dt[[d]]) & dt[[d]] == 1
    treat_var <- d
  } else {
    treated_mask <- rep(FALSE, nrow(dt))
    treat_var <- NULL
  }

  cluster_has_treated <- tapply(
    treated_mask,
    dt[[cluster_var]],
    any
  )
  n_treated_clusters <- as.integer(sum(cluster_has_treated))
  n_control_clusters <- as.integer(n_clusters - n_treated_clusters)

  if (is.null(treat_var)) {
    treatment_per_cluster <- rep.int(1L, n_clusters)
    names(treatment_per_cluster) <- cluster_ids
  } else {
    treatment_dt <- dt[
      ,
      .(n_unique = data.table::uniqueN(get(treat_var), na.rm = TRUE)),
      by = cluster_var
    ]
    treatment_per_cluster <- stats::setNames(
      as.integer(treatment_dt$n_unique),
      as.character(treatment_dt[[cluster_var]])
    )
  }

  treatment_varies_within <- any(treatment_per_cluster > 1L)
  n_clusters_with_variation <- as.integer(sum(treatment_per_cluster > 1L))

  units_per_cluster_dt <- dt[
    ,
    .(n_units = data.table::uniqueN(get(ivar))),
    by = cluster_var
  ]
  n_unique_units <- data.table::uniqueN(dt[[ivar]])
  is_nested_in_unit <- all(units_per_cluster_dt$n_units == 1L) &&
    n_clusters > n_unique_units
  level_relative_to_unit <- .determine_cluster_level(dt, ivar, cluster_var)

  if (length(cluster_sizes) <= 1L) {
    cv <- 0
  } else {
    cv <- stats::sd(cluster_sizes) / mean(cluster_sizes)
  }

  if (n_clusters > 0L) {
    balance_score <- min(
      min(n_treated_clusters, n_control_clusters) / (n_clusters / 2),
      1
    )
  } else {
    balance_score <- 0
  }
  cv_score <- max(0, 1 - cv / 2)
  reliability_score <- 0.5 * min(n_clusters / 50, 1) +
    0.3 * balance_score +
    0.2 * cv_score

  is_valid <- !is_nested_in_unit &&
    n_clusters >= 2L &&
    !identical(level_relative_to_unit, "lower")

  list(
    var_name = cluster_var,
    n_clusters = n_clusters,
    cluster_sizes = cluster_sizes,
    min_size = if (n_clusters > 0L) as.integer(min(cluster_sizes)) else 0L,
    max_size = if (n_clusters > 0L) as.integer(max(cluster_sizes)) else 0L,
    mean_size = if (n_clusters > 0L) mean(cluster_sizes) else 0,
    median_size = if (n_clusters > 0L) stats::median(cluster_sizes) else 0,
    cv = as.numeric(cv),
    n_treated_clusters = n_treated_clusters,
    n_control_clusters = n_control_clusters,
    treatment_varies_within = treatment_varies_within,
    is_nested_in_unit = is_nested_in_unit,
    level_relative_to_unit = level_relative_to_unit,
    units_per_cluster = mean(units_per_cluster_dt$n_units),
    n_clusters_with_variation = n_clusters_with_variation,
    balance_score = as.numeric(balance_score),
    cv_score = as.numeric(cv_score),
    reliability_score = as.numeric(reliability_score),
    is_valid = is_valid
  )
}


#' Detect the treatment variation level
#'
#' @param data data.frame or data.table.
#' @param ivar character(1). Unit identifier.
#' @param potential_cluster_vars character vector of candidate cluster variables.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param d character(1) or NULL. Treatment indicator.
#' @return character(1). First candidate level where treatment does not vary
#'   within groups; otherwise returns \code{ivar}.
#' @keywords internal
.detect_treatment_variation_level <- function(
    data, ivar, potential_cluster_vars, gvar, d
) {
  treat_var <- if (!is.null(gvar)) gvar else d
  if (is.null(treat_var)) {
    return(ivar)
  }

  dt <- data.table::as.data.table(data)
  n_levels <- vapply(
    potential_cluster_vars,
    function(var) data.table::uniqueN(dt[[var]], na.rm = TRUE),
    integer(1)
  )
  sorted_vars <- potential_cluster_vars[order(n_levels)]

  for (cluster_var in sorted_vars) {
    treatment_per_group <- dt[
      ,
      .(n_unique = data.table::uniqueN(get(treat_var), na.rm = TRUE)),
      by = cluster_var
    ]

    if (all(treatment_per_group$n_unique == 1L)) {
      return(cluster_var)
    }
  }

  ivar
}


#' Generate an internal clustering recommendation
#'
#' @param cluster_structure named list of cluster diagnostics.
#' @param treatment_level character(1). Detected treatment variation level.
#' @return list with \code{recommended_var}, \code{reason}, and \code{warnings}.
#' @keywords internal
.generate_clustering_recommendation <- function(
    cluster_structure, treatment_level
) {
  warnings <- character(0)
  valid_options <- Filter(function(x) isTRUE(x$is_valid), cluster_structure)

  if (length(valid_options) == 0L) {
    return(list(
      recommended_var = NULL,
      reason = "No valid clustering options available.",
      warnings = c(
        "All potential cluster variables are invalid (nested within units or < 2 clusters)."
      )
    ))
  }

  if (treatment_level %in% names(valid_options)) {
    treatment_stats <- valid_options[[treatment_level]]
    if (isTRUE(treatment_stats$n_clusters >= 20L)) {
      return(list(
        recommended_var = treatment_level,
        reason = paste0(
          "Clustering at treatment variation level (",
          treatment_level,
          ") with ",
          treatment_stats$n_clusters,
          " clusters."
        ),
        warnings = warnings
      ))
    }

    warnings <- c(
      warnings,
      paste0(
        "Treatment varies at ",
        treatment_level,
        " level but only ",
        treatment_stats$n_clusters,
        " clusters available."
      )
    )
  }

  option_names <- names(valid_options)
  option_scores <- vapply(
    valid_options,
    function(x) x$reliability_score,
    numeric(1)
  )
  option_large_n <- vapply(
    valid_options,
    function(x) as.integer(x$n_clusters >= 20L),
    integer(1)
  )
  ranked <- option_names[order(option_large_n, option_scores, decreasing = TRUE)]
  best_var <- ranked[[1L]]
  best_stats <- valid_options[[best_var]]

  if (best_stats$n_clusters < 20L) {
    warnings <- c(
      warnings,
      paste0(
        "Recommended clustering has only ",
        best_stats$n_clusters,
        " clusters. Consider wild cluster bootstrap for reliable inference."
      )
    )
  }

  if (isTRUE(best_stats$treatment_varies_within)) {
    warnings <- c(
      warnings,
      paste0(
        "Treatment varies within ",
        best_stats$n_clusters_with_variation,
        " clusters. Standard errors may be conservative."
      )
    )
  }

  list(
    recommended_var = best_var,
    reason = paste0(
      "Best available option with ",
      best_stats$n_clusters,
      " clusters (reliability score: ",
      formatC(best_stats$reliability_score, format = "f", digits = 2L),
      ")."
    ),
    warnings = warnings
  )
}


.generate_recommendation_reasons <- function(var, stats, diag) {
  reasons <- character(0)

  if (identical(var, diag$treatment_variation_level)) {
    reasons <- c(
      reasons,
      paste0(
        "Treatment varies at ",
        var,
        " level - clustering at this level is appropriate"
      )
    )
  }

  if (stats$n_clusters >= 30L) {
    reasons <- c(
      reasons,
      paste0(
        "Sufficient clusters (",
        stats$n_clusters,
        ") for reliable inference"
      )
    )
  } else if (stats$n_clusters >= 20L) {
    reasons <- c(
      reasons,
      paste0(
        "Adequate clusters (",
        stats$n_clusters,
        ") for inference"
      )
    )
  } else {
    reasons <- c(
      reasons,
      paste0(
        "Limited clusters (",
        stats$n_clusters,
        ") - consider wild bootstrap"
      )
    )
  }

  if (stats$n_treated_clusters > 0L && stats$n_control_clusters > 0L) {
    balance_ratio <- min(stats$n_treated_clusters, stats$n_control_clusters) /
      max(stats$n_treated_clusters, stats$n_control_clusters)
    if (balance_ratio > 0.5) {
      reasons <- c(
        reasons,
        paste0(
          "Good balance between treated (",
          stats$n_treated_clusters,
          ") and control (",
          stats$n_control_clusters,
          ") clusters"
        )
      )
    }
  }

  if (identical(stats$level_relative_to_unit, "higher")) {
    reasons <- c(
      reasons,
      "Clustering at higher level than unit of observation"
    )
  }

  reasons
}


.get_alternative_reason <- function(stats) {
  if (stats$n_clusters < 10L) {
    return("too few clusters")
  }
  if (stats$n_clusters < 20L) {
    return("marginal cluster count")
  }
  if (isTRUE(stats$treatment_varies_within)) {
    return("treatment varies within clusters")
  }
  paste0(
    "reliability score: ",
    formatC(stats$reliability_score, format = "f", digits = 2L)
  )
}


.generate_clustering_warnings <- function(stats) {
  warnings <- character(0)

  if (stats$n_clusters < 10L) {
    warnings <- c(
      warnings,
      paste0(
        "Only ",
        stats$n_clusters,
        " clusters - cluster-robust inference may be unreliable"
      )
    )
  } else if (stats$n_clusters < 20L) {
    warnings <- c(
      warnings,
      paste0(
        "Only ",
        stats$n_clusters,
        " clusters - consider wild cluster bootstrap"
      )
    )
  }

  if (isTRUE(stats$cv > 1.0)) {
    warnings <- c(
      warnings,
      paste0(
        "Highly variable cluster sizes (CV=",
        formatC(stats$cv, format = "f", digits = 2L),
        ")"
      )
    )
  }

  if (isTRUE(stats$treatment_varies_within)) {
    warnings <- c(
      warnings,
      paste0(
        "Treatment varies within ",
        stats$n_clusters_with_variation,
        " clusters"
      )
    )
  }

  warnings
}


.emit_clustering_diagnosis_summary <- function(x) {
  cat("Clustering Diagnostics\n")
  cat(
    sprintf(
      "Treatment variation level: %s\n",
      x$treatment_variation_level
    )
  )
  cat(
    sprintf(
      "Recommended cluster variable: %s\n",
      if (is.null(x$recommended_cluster_var)) "None" else x$recommended_cluster_var
    )
  )
  cat(
    sprintf(
      "Recommendation reason: %s\n",
      x$recommendation_reason
    )
  )

  if (length(x$warnings) > 0L) {
    cat("Warnings:\n")
    for (warning_text in x$warnings) {
      cat(sprintf("  - %s\n", warning_text))
    }
  }

  invisible(x)
}


.emit_clustering_recommendation_summary <- function(x) {
  cat("Clustering Recommendation\n")
  cat(
    sprintf(
      "Recommended cluster variable: %s\n",
      if (is.null(x$recommended_var)) "None" else x$recommended_var
    )
  )

  if (!is.null(x$recommended_var)) {
    cat(
      sprintf(
        "Clusters: %d | Treated: %d | Control: %d\n",
        x$n_clusters,
        x$n_treated_clusters,
        x$n_control_clusters
      )
    )
    cat(
      sprintf(
        "Confidence: %s\n",
        formatC(x$confidence, format = "f", digits = 2L)
      )
    )
  }

  if (length(x$reasons) > 0L) {
    cat("Reasons:\n")
    for (reason_text in x$reasons) {
      cat(sprintf("  - %s\n", reason_text))
    }
  }

  if (length(x$alternatives) > 0L) {
    cat("Alternatives:\n")
    for (alternative in x$alternatives) {
      cat(
        sprintf(
          "  - %s (%s)\n",
          alternative$var,
          alternative$reason
        )
      )
    }
  }

  if (length(x$warnings) > 0L) {
    cat("Warnings:\n")
    for (warning_text in x$warnings) {
      cat(sprintf("  - %s\n", warning_text))
    }
  }

  if (isTRUE(x$use_wild_bootstrap) && !is.null(x$wild_bootstrap_reason)) {
    cat(
      sprintf(
        "Wild bootstrap: %s\n",
        x$wild_bootstrap_reason
      )
    )
  }

  invisible(x)
}


.emit_clustering_consistency_summary <- function(x) {
  cat("Clustering Consistency Check\n")
  cat(
    sprintf(
      "Consistent: %s\n",
      if (isTRUE(x$is_consistent)) "YES" else "NO"
    )
  )
  cat(
    sprintf(
      "Treatment variation level: %s\n",
      x$treatment_variation_level
    )
  )
  cat(
    sprintf(
      "Cluster level: %s\n",
      x$cluster_level
    )
  )
  cat(
    sprintf(
      "Inconsistent clusters: %d of %d (%.1f%%)\n",
      x$n_inconsistent,
      x$n_clusters,
      x$pct_inconsistent
    )
  )
  cat(sprintf("Recommendation: %s\n", x$recommendation))
  invisible(x)
}


#' Diagnose clustering structure for inference choices
#'
#' @param data data.frame or data.table in long format.
#' @param ivar character(1). Unit identifier.
#' @param potential_cluster_vars character vector of candidate cluster variables.
#' @param gvar character(1) or NULL. Cohort variable for staggered designs.
#' @param d character(1) or NULL. Treatment indicator for common-timing designs.
#' @param verbose logical(1). If \code{TRUE}, print a concise diagnosis summary.
#' @return An object of class \code{lwdid_clustering_diagnosis}.
#' @export
diagnose_clustering <- function(
    data, ivar, potential_cluster_vars, gvar = NULL, d = NULL, verbose = TRUE
) {
  .validate_clustering_inputs(data, ivar, potential_cluster_vars, gvar, d)

  cluster_structure <- lapply(
    potential_cluster_vars,
    function(cluster_var) .analyze_cluster_var(data, ivar, cluster_var, gvar, d)
  )
  names(cluster_structure) <- potential_cluster_vars

  treatment_level <- .detect_treatment_variation_level(
    data,
    ivar,
    potential_cluster_vars,
    gvar,
    d
  )
  recommendation <- .generate_clustering_recommendation(
    cluster_structure,
    treatment_level
  )

  result <- structure(
    list(
      cluster_structure = cluster_structure,
      recommended_cluster_var = recommendation$recommended_var,
      recommendation_reason = recommendation$reason,
      treatment_variation_level = treatment_level,
      warnings = recommendation$warnings
    ),
    class = "lwdid_clustering_diagnosis"
  )

  if (isTRUE(verbose)) {
    .emit_clustering_diagnosis_summary(result)
  }

  result
}


#' Print method for clustering diagnosis objects
#'
#' @param x An object of class \code{lwdid_clustering_diagnosis}.
#' @param ... Additional arguments (ignored).
#' @return \code{x} invisibly.
#' @export
print.lwdid_clustering_diagnosis <- function(x, ...) {
  .emit_clustering_diagnosis_summary(x)
}


#' Recommend a clustering level with detailed guidance
#'
#' @param data data.frame or data.table in long format.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time identifier.
#' @param potential_cluster_vars character vector of candidate cluster variables.
#' @param gvar character(1) or NULL. Cohort variable for staggered designs.
#' @param d character(1) or NULL. Treatment indicator for common-timing designs.
#' @param min_clusters integer(1). Threshold for wild cluster bootstrap advice.
#' @param verbose logical(1). If \code{TRUE}, print a concise recommendation.
#' @return An object of class \code{lwdid_clustering_recommendation}.
#' @export
recommend_clustering <- function(
    data, ivar, tvar, potential_cluster_vars, gvar = NULL, d = NULL,
    min_clusters = 20L, verbose = TRUE
) {
  .validate_clustering_inputs(data, ivar, potential_cluster_vars, gvar, d)

  if (!is.character(tvar) || length(tvar) != 1L || !tvar %in% names(data)) {
    stop(sprintf("Time variable '%s' not found in data", tvar), call. = FALSE)
  }

  diag <- diagnose_clustering(
    data = data,
    ivar = ivar,
    potential_cluster_vars = potential_cluster_vars,
    gvar = gvar,
    d = d,
    verbose = FALSE
  )

  valid_options <- Filter(function(x) isTRUE(x$is_valid), diag$cluster_structure)
  if (length(valid_options) == 0L || is.null(diag$recommended_cluster_var)) {
    result <- structure(
      list(
        recommended_var = NULL,
        n_clusters = NA_integer_,
        n_treated_clusters = NA_integer_,
        n_control_clusters = NA_integer_,
        confidence = 0,
        reasons = "No valid clustering options available.",
        alternatives = list(),
        warnings = diag$warnings,
        use_wild_bootstrap = FALSE,
        wild_bootstrap_reason = NULL
      ),
      class = "lwdid_clustering_recommendation"
    )

    if (isTRUE(verbose)) {
      .emit_clustering_recommendation_summary(result)
    }

    return(result)
  }

  recommended_var <- diag$recommended_cluster_var
  best_stats <- diag$cluster_structure[[recommended_var]]
  reasons <- .generate_recommendation_reasons(recommended_var, best_stats, diag)

  option_names <- names(valid_options)
  option_scores <- vapply(
    valid_options,
    function(x) x$reliability_score,
    numeric(1)
  )
  option_large_n <- vapply(
    valid_options,
    function(x) as.integer(x$n_clusters >= 20L),
    integer(1)
  )
  ranked <- option_names[order(option_large_n, option_scores, decreasing = TRUE)]
  alternative_names <- setdiff(ranked, recommended_var)
  alternatives <- lapply(
    utils::head(alternative_names, 2L),
    function(var) {
      stats <- diag$cluster_structure[[var]]
      list(
        var = var,
        n_clusters = stats$n_clusters,
        reliability_score = stats$reliability_score,
        reason = .get_alternative_reason(stats)
      )
    }
  )

  use_wild_bootstrap <- isTRUE(best_stats$n_clusters < min_clusters)
  wild_bootstrap_reason <- NULL
  if (use_wild_bootstrap) {
    wild_bootstrap_reason <- paste0(
      "Only ",
      best_stats$n_clusters,
      " clusters available (< ",
      min_clusters,
      "). Wild cluster bootstrap recommended for reliable inference."
    )
  }

  result <- structure(
    list(
      recommended_var = recommended_var,
      n_clusters = as.integer(best_stats$n_clusters),
      n_treated_clusters = as.integer(best_stats$n_treated_clusters),
      n_control_clusters = as.integer(best_stats$n_control_clusters),
      confidence = as.numeric(best_stats$reliability_score),
      reasons = reasons,
      alternatives = alternatives,
      warnings = .generate_clustering_warnings(best_stats),
      use_wild_bootstrap = use_wild_bootstrap,
      wild_bootstrap_reason = wild_bootstrap_reason
    ),
    class = "lwdid_clustering_recommendation"
  )

  if (isTRUE(verbose)) {
    .emit_clustering_recommendation_summary(result)
  }

  result
}


#' Print method for clustering recommendation objects
#'
#' @param x An object of class \code{lwdid_clustering_recommendation}.
#' @param ... Additional arguments (ignored).
#' @return \code{x} invisibly.
#' @export
print.lwdid_clustering_recommendation <- function(x, ...) {
  .emit_clustering_recommendation_summary(x)
}


#' Check consistency between treatment assignment and clustering choice
#'
#' @param data data.frame or data.table in long format.
#' @param ivar character(1). Unit identifier.
#' @param cluster_var character(1). Candidate clustering variable to validate.
#' @param gvar character(1) or NULL. Cohort variable for staggered designs.
#' @param d character(1) or NULL. Treatment indicator for common-timing designs.
#' @param verbose logical(1). If \code{TRUE}, print a concise consistency report.
#' @return A named list describing clustering consistency.
#' @export
check_clustering_consistency <- function(
    data, ivar, cluster_var, gvar = NULL, d = NULL, verbose = TRUE
) {
  .validate_clustering_inputs(data, ivar, cluster_var, gvar, d)

  dt <- data.table::as.data.table(data)
  treat_var <- if (!is.null(gvar)) gvar else d
  cluster_treatment <- dt[
    ,
    .(n_unique = data.table::uniqueN(get(treat_var), na.rm = TRUE)),
    by = cluster_var
  ]
  inconsistent_clusters <- as.character(
    cluster_treatment[[cluster_var]][cluster_treatment$n_unique > 1L]
  )
  n_clusters <- nrow(cluster_treatment)
  n_inconsistent <- length(inconsistent_clusters)
  pct_inconsistent <- if (n_clusters > 0L) {
    n_inconsistent / n_clusters * 100
  } else {
    0
  }

  treatment_variation_level <- .detect_treatment_variation_level(
    data,
    ivar,
    c(cluster_var),
    gvar,
    d
  )
  cluster_level <- .determine_cluster_level(data, ivar, cluster_var)
  is_consistent <- isTRUE(pct_inconsistent < 5) &&
    cluster_level %in% c("same", "higher")

  if (is_consistent) {
    recommendation <- "Clustering choice is appropriate."
  } else if (pct_inconsistent >= 5) {
    recommendation <- paste0(
      "Treatment varies within ",
      formatC(pct_inconsistent, format = "f", digits = 1L),
      "% of clusters. Consider clustering at a higher level where treatment is constant."
    )
  } else {
    recommendation <- paste0(
      "Cluster level (",
      cluster_level,
      ") may be inappropriate. Consider clustering at the treatment variation level (",
      treatment_variation_level,
      ")."
    )
  }

  details <- paste0(
    "Analyzed ",
    n_clusters,
    " clusters.\n",
    "Treatment varies within ",
    n_inconsistent,
    " clusters (",
    formatC(pct_inconsistent, format = "f", digits = 1L),
    "%).\n",
    "Treatment variation level: ",
    treatment_variation_level,
    "\n",
    "Cluster level: ",
    cluster_level
  )

  result <- list(
    is_consistent = is_consistent,
    treatment_variation_level = treatment_variation_level,
    cluster_level = cluster_level,
    n_clusters = as.integer(n_clusters),
    n_inconsistent = as.integer(n_inconsistent),
    pct_inconsistent = as.numeric(pct_inconsistent),
    inconsistent_clusters = inconsistent_clusters,
    recommendation = recommendation,
    details = details
  )

  if (isTRUE(verbose)) {
    .emit_clustering_consistency_summary(result)
  }

  result
}


.parse_normal_distribution <- function(spec, default_mean, default_sd) {
  if (is.null(spec) || !is.character(spec) || length(spec) != 1L) {
    return(list(mean = default_mean, sd = default_sd))
  }

  matches <- regexec(
    "^normal\\(([-0-9.]+),\\s*([-0-9.]+)\\)$",
    spec
  )
  parsed <- regmatches(spec, matches)[[1L]]

  if (length(parsed) != 3L) {
    return(list(mean = default_mean, sd = default_sd))
  }

  list(
    mean = as.numeric(parsed[[2L]]),
    sd = as.numeric(parsed[[3L]])
  )
}


.draw_cluster_sizes <- function(scenario) {
  if (!is.null(scenario[["obs_per_cluster"]])) {
    return(rep.int(as.integer(scenario[["obs_per_cluster"]]), as.integer(scenario[["G"]])))
  }

  if (identical(scenario[["cluster_size_draw"]], "randint(10, 101)")) {
    return(sample(10:100, size = as.integer(scenario[["G"]]), replace = TRUE))
  }

  stop(
    sprintf(
      "Unsupported clustering Monte Carlo size design for scenario '%s'.",
      scenario[["id"]]
    ),
    call. = FALSE
  )
}


.draw_cluster_treatment <- function(g, probability) {
  repeat {
    treatment <- stats::rbinom(as.integer(g), size = 1L, prob = probability)
    if (any(treatment == 0L) && any(treatment == 1L)) {
      return(treatment)
    }
  }
}


.simulate_clustering_monte_carlo_panel <- function(scenario, shared_dgp = NULL) {
  intercept <- if (!is.null(shared_dgp[["intercept"]])) {
    as.numeric(shared_dgp[["intercept"]])
  } else {
    10
  }
  cluster_dist <- .parse_normal_distribution(
    shared_dgp[["cluster_effect_distribution"]],
    default_mean = 0,
    default_sd = 2
  )
  error_dist <- .parse_normal_distribution(
    shared_dgp[["idiosyncratic_error_distribution"]],
    default_mean = 0,
    default_sd = 1
  )
  treatment_probability <- if (!is.null(scenario[["treatment_probability"]])) {
    as.numeric(scenario[["treatment_probability"]])
  } else if (!is.null(shared_dgp[["treatment_assignment"]][["default_probability"]])) {
    as.numeric(shared_dgp[["treatment_assignment"]][["default_probability"]])
  } else {
    0.5
  }

  cluster_sizes <- .draw_cluster_sizes(scenario)
  g <- as.integer(scenario[["G"]])
  cluster_ids <- rep.int(seq_len(g), times = cluster_sizes)
  cluster_treatment <- .draw_cluster_treatment(g, treatment_probability)
  treatment <- rep.int(cluster_treatment, times = cluster_sizes)
  cluster_effects <- rep.int(
    stats::rnorm(g, mean = cluster_dist$mean, sd = cluster_dist$sd),
    times = cluster_sizes
  )
  errors <- stats::rnorm(
    length(cluster_ids),
    mean = error_dist$mean,
    sd = error_dist$sd
  )
  outcome <- intercept +
    as.numeric(scenario[["true_tau"]]) * treatment +
    cluster_effects +
    errors

  data.frame(
    Y = outcome,
    D = treatment,
    cluster = cluster_ids
  )
}


.prepare_clustering_monte_carlo_precomp <- function(d, cluster) {
  cluster_index <- as.integer(as.factor(cluster))
  x <- cbind("(Intercept)" = 1, D = as.numeric(d))
  xtx_inv <- solve(crossprod(x))
  n_obs <- nrow(x)
  g <- length(unique(cluster_index))
  correction <- (g / (g - 1)) * ((n_obs - 1) / (n_obs - ncol(x)))

  list(
    X = x,
    XtX_inv = xtx_inv,
    cluster_index = cluster_index,
    cluster_rows = split(seq_len(n_obs), cluster_index),
    G = g,
    correction = correction,
    bread_row_d = as.numeric(xtx_inv[2L, ])
  )
}


.fit_clustering_monte_carlo_model <- function(y, precomp) {
  beta <- drop(precomp$XtX_inv %*% crossprod(precomp$X, y))
  fitted <- drop(precomp$X %*% beta)
  residuals <- y - fitted

  cluster_scores <- vapply(
    precomp$cluster_rows,
    function(idx) {
      drop(crossprod(precomp$X[idx, , drop = FALSE], residuals[idx]))
    },
    numeric(2)
  )

  if (is.null(dim(cluster_scores))) {
    cluster_scores <- matrix(cluster_scores, nrow = 2L)
  }

  meat <- cluster_scores %*% t(cluster_scores)
  var_d <- precomp$correction * as.numeric(
    t(precomp$bread_row_d) %*% meat %*% precomp$bread_row_d
  )
  se <- sqrt(max(var_d, 0))

  list(
    att = as.numeric(beta[[2L]]),
    se = se,
    fitted = fitted,
    residuals = residuals
  )
}


.fit_clustering_monte_carlo_model_batch <- function(y_matrix, precomp) {
  y_matrix <- as.matrix(y_matrix)
  beta_mat <- precomp$XtX_inv %*% crossprod(precomp$X, y_matrix)
  residual_mat <- y_matrix - precomp$X %*% beta_mat

  n_boot <- ncol(y_matrix)
  meat11 <- numeric(n_boot)
  meat12 <- numeric(n_boot)
  meat22 <- numeric(n_boot)

  for (idx in precomp$cluster_rows) {
    cluster_scores <- crossprod(
      precomp$X[idx, , drop = FALSE],
      residual_mat[idx, , drop = FALSE]
    )
    meat11 <- meat11 + cluster_scores[1L, ]^2
    meat12 <- meat12 + cluster_scores[1L, ] * cluster_scores[2L, ]
    meat22 <- meat22 + cluster_scores[2L, ]^2
  }

  a1 <- precomp$bread_row_d[[1L]]
  a2 <- precomp$bread_row_d[[2L]]
  var_d <- precomp$correction * (
    a1^2 * meat11 +
      2 * a1 * a2 * meat12 +
      a2^2 * meat22
  )

  list(
    att = as.numeric(beta_mat[2L, ]),
    se = sqrt(pmax(var_d, 0))
  )
}


.generate_all_rademacher_weights <- function(n_clusters) {
  weights <- expand.grid(
    rep(list(c(-1, 1)), n_clusters),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  data.matrix(weights)
}


.resolve_wild_bootstrap_weights <- function(
    n_clusters, requested_n_bootstrap, weight_type, full_enumeration = NULL
) {
  supported_weight_types <- c("rademacher", "mammen", "webb")
  if (!weight_type %in% supported_weight_types) {
    stop(
      sprintf(
        "Unsupported wild bootstrap weight_type '%s'. Must be one of: %s.",
        weight_type,
        paste(sprintf("'%s'", supported_weight_types), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (is.null(full_enumeration)) {
    full_enumeration <- n_clusters <= 12L && identical(weight_type, "rademacher")
  }

  if (isTRUE(full_enumeration)) {
    if (!identical(weight_type, "rademacher")) {
      stop(
        "Full enumeration is currently only available for weight_type = 'rademacher'.",
        call. = FALSE
      )
    }
    weights <- .generate_all_rademacher_weights(n_clusters)
  } else {
    weights <- t(vapply(
      seq_len(as.integer(requested_n_bootstrap)),
      function(i) {
        .generate_weights(n_clusters, weight_type)
      },
      numeric(as.integer(n_clusters))
    ))
  }

  list(
    weights = weights,
    requested_n_bootstrap = as.integer(requested_n_bootstrap),
    actual_n_bootstrap = nrow(weights),
    full_enumeration = isTRUE(full_enumeration)
  )
}


.run_wild_cluster_bootstrap_ci <- function(
    data, requested_n_bootstrap, weight_type, alpha, impose_null,
    full_enumeration = NULL
) {
  precomp <- .prepare_clustering_monte_carlo_precomp(data$D, data$cluster)
  original_fit <- .fit_clustering_monte_carlo_model(data$Y, precomp)

  if (!is.finite(original_fit$se) || original_fit$se <= 0) {
    return(list(
      att = original_fit$att,
      se = original_fit$se,
      t_stat_original = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      pvalue = NA_real_,
      requested_n_bootstrap = as.integer(requested_n_bootstrap),
      actual_n_bootstrap = as.integer(requested_n_bootstrap),
      weight_type = weight_type,
      impose_null = isTRUE(impose_null),
      full_enumeration = FALSE,
      ci_method = if (isTRUE(impose_null)) "percentile_t" else "percentile"
    ))
  }

  if (isTRUE(impose_null)) {
    fitted_base <- rep.int(mean(data$Y), nrow(data))
    residual_base <- data$Y - fitted_base
  } else {
    fitted_base <- original_fit$fitted
    residual_base <- original_fit$residuals
  }

  weight_info <- .resolve_wild_bootstrap_weights(
    n_clusters = precomp$G,
    requested_n_bootstrap = requested_n_bootstrap,
    weight_type = weight_type,
    full_enumeration = full_enumeration
  )

  obs_weights <- t(weight_info$weights[, precomp$cluster_index, drop = FALSE])
  y_boot <- matrix(
    fitted_base,
    nrow = nrow(data),
    ncol = nrow(weight_info$weights)
  ) + residual_base * obs_weights

  boot_fit <- .fit_clustering_monte_carlo_model_batch(y_boot, precomp)
  valid <- is.finite(boot_fit$att) & is.finite(boot_fit$se) & boot_fit$se > 0

  if (!any(valid)) {
    return(list(
      att = original_fit$att,
      se = original_fit$se,
      t_stat_original = original_fit$att / original_fit$se,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      pvalue = NA_real_,
      requested_n_bootstrap = weight_info$requested_n_bootstrap,
      actual_n_bootstrap = weight_info$actual_n_bootstrap,
      weight_type = weight_type,
      impose_null = isTRUE(impose_null),
      full_enumeration = weight_info$full_enumeration,
      ci_method = if (isTRUE(impose_null)) "percentile_t" else "percentile"
    ))
  }

  att_valid <- boot_fit$att[valid]
  t_stats <- boot_fit$att[valid] / boot_fit$se[valid]
  t_stat_original <- original_fit$att / original_fit$se
  if (isTRUE(impose_null)) {
    t_abs_crit <- as.numeric(
      stats::quantile(abs(t_stats), probs = 1 - alpha, names = FALSE, type = 7)
    )
    ci_lower <- original_fit$att - t_abs_crit * original_fit$se
    ci_upper <- original_fit$att + t_abs_crit * original_fit$se
  } else {
    ci_lower <- as.numeric(
      stats::quantile(att_valid, probs = alpha / 2, names = FALSE, type = 7)
    )
    ci_upper <- as.numeric(
      stats::quantile(att_valid, probs = 1 - alpha / 2, names = FALSE, type = 7)
    )
  }

  list(
    att = original_fit$att,
    se = original_fit$se,
    t_stat_original = t_stat_original,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    pvalue = mean(abs(t_stats) >= abs(t_stat_original)),
    requested_n_bootstrap = weight_info$requested_n_bootstrap,
    actual_n_bootstrap = weight_info$actual_n_bootstrap,
    weight_type = weight_type,
    impose_null = isTRUE(impose_null),
    full_enumeration = weight_info$full_enumeration,
    ci_method = if (isTRUE(impose_null)) "percentile_t" else "percentile"
  )
}


.evaluate_clustering_monte_carlo_standard_ci <- function(data, true_tau, alpha) {
  fit <- stats::lm(Y ~ D, data = data)
  coefficients <- stats::coef(fit)

  if (!("D" %in% names(coefficients)) || !is.finite(coefficients[["D"]])) {
    return(NULL)
  }

  vcov_result <- withCallingHandlers(
    compute_cluster_vce(fit, cluster = data$cluster, type = "HC1"),
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  )
  se <- sqrt(diag(vcov_result$vcov))[["D"]]
  if (!is.finite(se) || se <= 0) {
    return(NULL)
  }

  t_crit <- stats::qt(1 - alpha / 2, df = vcov_result$df)
  ci_lower <- coefficients[["D"]] - t_crit * se
  ci_upper <- coefficients[["D"]] + t_crit * se

  list(
    att = coefficients[["D"]],
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    covers = isTRUE(ci_lower <= true_tau && true_tau <= ci_upper)
  )
}


.run_clustering_monte_carlo_wild_bootstrap_scenario <- function(
    scenario, shared_dgp = NULL
) {
  set.seed(as.integer(scenario[["seed"]]))

  n_target <- as.integer(scenario[["n_simulations"]])
  alpha <- if (!is.null(scenario[["alpha"]])) {
    as.numeric(scenario[["alpha"]])
  } else {
    0.05
  }

  requested_n_bootstrap <- if (!is.null(scenario[["requested_n_bootstrap"]])) {
    as.integer(scenario[["requested_n_bootstrap"]])
  } else {
    as.integer(scenario[["n_bootstrap"]])
  }

  metric_hits <- logical(n_target)
  replication_rows <- vector("list", n_target)
  valid_replications <- 0L
  attempts <- 0L
  max_attempts <- n_target * 20L
  last_result <- NULL

  while (valid_replications < n_target) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop(
        sprintf(
          "Clustering Monte Carlo scenario '%s' exceeded %d attempts.",
          scenario[["id"]],
          max_attempts
        ),
        call. = FALSE
      )
    }

    sim_data <- .simulate_clustering_monte_carlo_panel(scenario, shared_dgp)
    wild_result <- .run_wild_cluster_bootstrap_ci(
      data = sim_data,
      requested_n_bootstrap = requested_n_bootstrap,
      weight_type = scenario[["weight_type"]],
      alpha = alpha,
      impose_null = scenario[["impose_null"]]
    )

    if (!is.finite(wild_result$ci_lower) || !is.finite(wild_result$ci_upper)) {
      next
    }

    valid_replications <- valid_replications + 1L
    metric_hits[[valid_replications]] <- isTRUE(
      wild_result$ci_lower <= as.numeric(scenario[["true_tau"]]) &&
        as.numeric(scenario[["true_tau"]]) <= wild_result$ci_upper
    )
    last_result <- wild_result
  }

  list(
    scenario_id = scenario[["id"]],
    n_simulations = n_target,
    n_attempts = attempts,
    metric_name = "coverage_rate",
    metric_value = mean(metric_hits),
    requested_n_bootstrap = last_result$requested_n_bootstrap,
    actual_n_bootstrap = last_result$actual_n_bootstrap,
    weight_type = last_result$weight_type,
    impose_null = last_result$impose_null,
    full_enumeration = last_result$full_enumeration,
    ci_method = last_result$ci_method
  )
}


.run_clustering_monte_carlo_wild_bootstrap_decision_scenario <- function(
    scenario, shared_dgp = NULL
) {
  set.seed(as.integer(scenario[["seed"]]))

  n_target <- as.integer(scenario[["n_simulations"]])
  alpha <- if (!is.null(scenario[["alpha"]])) {
    as.numeric(scenario[["alpha"]])
  } else {
    0.05
  }
  requested_n_bootstrap <- if (!is.null(scenario[["requested_n_bootstrap"]])) {
    as.integer(scenario[["requested_n_bootstrap"]])
  } else {
    as.integer(scenario[["n_bootstrap"]])
  }
  pvalue_threshold <- if (!is.null(scenario[["pvalue_threshold"]])) {
    as.numeric(scenario[["pvalue_threshold"]])
  } else {
    alpha
  }
  decision_direction <- scenario[["decision_direction"]]
  if (!decision_direction %in% c("gt", "lt")) {
    stop(
      sprintf(
        "Unsupported decision_direction '%s' for clustering Monte Carlo scenario '%s'.",
        decision_direction,
        scenario[["id"]]
      ),
      call. = FALSE
    )
  }

  metric_hits <- logical(n_target)
  replication_rows <- vector("list", n_target)
  valid_replications <- 0L
  attempts <- 0L
  max_attempts <- n_target * 20L
  last_result <- NULL

  while (valid_replications < n_target) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop(
        sprintf(
          "Clustering Monte Carlo scenario '%s' exceeded %d attempts.",
          scenario[["id"]],
          max_attempts
        ),
        call. = FALSE
      )
    }

    sim_data <- .simulate_clustering_monte_carlo_panel(scenario, shared_dgp)
    wild_result <- .run_wild_cluster_bootstrap_ci(
      data = sim_data,
      requested_n_bootstrap = requested_n_bootstrap,
      weight_type = scenario[["weight_type"]],
      alpha = alpha,
      impose_null = scenario[["impose_null"]]
    )

    if (!is.finite(wild_result$pvalue)) {
      next
    }

    standard_result <- .evaluate_clustering_monte_carlo_standard_ci(
      data = sim_data,
      true_tau = as.numeric(scenario[["true_tau"]]),
      alpha = alpha
    )
    standard_t_pvalue <- NA_real_
    standard_t_reject <- NA
    if (!is.null(standard_result) &&
        is.finite(standard_result$att) &&
        is.finite(standard_result$se) &&
        standard_result$se > 0) {
      standard_t_stat <- standard_result$att / standard_result$se
      standard_t_pvalue <- 2 * stats::pt(
        -abs(standard_t_stat),
        df = as.integer(scenario[["G"]]) - 1L
      )
      standard_t_reject <- isTRUE(standard_t_pvalue < pvalue_threshold)
    }

    valid_replications <- valid_replications + 1L
    metric_hits[[valid_replications]] <- if (identical(decision_direction, "gt")) {
      wild_result$pvalue > pvalue_threshold
    } else {
      wild_result$pvalue < pvalue_threshold
    }
    replication_rows[[valid_replications]] <- list(
      attempt_id = as.integer(attempts),
      replication_id = as.integer(valid_replications),
      treated_cluster_count = as.integer(length(unique(sim_data$cluster[sim_data$D == 1]))),
      reject = metric_hits[[valid_replications]],
      pvalue = as.numeric(wild_result$pvalue),
      att = as.numeric(wild_result$att),
      original_se = as.numeric(wild_result$se),
      t_stat_original = as.numeric(wild_result$t_stat_original),
      standard_t_pvalue = as.numeric(standard_t_pvalue),
      standard_t_reject = standard_t_reject
    )
    last_result <- wild_result
  }

  treated_counts <- vapply(
    replication_rows,
    `[[`,
    integer(1L),
    "treated_cluster_count"
  )
  treated_count_levels <- sort(unique(treated_counts))
  treated_cluster_profile <- lapply(
    treated_count_levels,
    function(count) {
      match_idx <- treated_counts == count
      bucket_rows <- replication_rows[match_idx]
      rejection_count <- sum(vapply(
        bucket_rows,
        `[[`,
        logical(1L),
        "reject"
      ))
      n_replications <- sum(match_idx)
      bucket_pvalues <- vapply(
        bucket_rows,
        `[[`,
        numeric(1L),
        "pvalue"
      )
      bucket_abs_att <- abs(vapply(
        bucket_rows,
        `[[`,
        numeric(1L),
        "att"
      ))
      bucket_original_se <- vapply(
        bucket_rows,
        `[[`,
        numeric(1L),
        "original_se"
      )
      bucket_att_to_se_ratio <- bucket_abs_att / bucket_original_se
      bucket_abs_t_stats <- abs(vapply(
        bucket_rows,
        `[[`,
        numeric(1L),
        "t_stat_original"
      ))

      list(
        treated_cluster_count = as.integer(count),
        n_replications = as.integer(n_replications),
        rejection_count = as.integer(rejection_count),
        rejection_rate = rejection_count / n_replications,
        avg_abs_att = mean(bucket_abs_att),
        median_abs_att = stats::median(bucket_abs_att),
        avg_original_se = mean(bucket_original_se),
        median_original_se = stats::median(bucket_original_se),
        avg_abs_att_to_se_ratio = mean(bucket_att_to_se_ratio),
        median_abs_att_to_se_ratio = stats::median(bucket_att_to_se_ratio),
        min_abs_att_to_se_ratio = min(bucket_att_to_se_ratio),
        max_abs_att_to_se_ratio = max(bucket_att_to_se_ratio),
        avg_pvalue = mean(bucket_pvalues),
        median_pvalue = stats::median(bucket_pvalues),
        avg_abs_t_stat = mean(bucket_abs_t_stats),
        max_abs_t_stat = max(bucket_abs_t_stats),
        near_threshold_count = as.integer(sum(
          bucket_pvalues >= 0.04 & bucket_pvalues <= 0.06
        ))
      )
    }
  )

  metric_name <- if (!is.null(scenario[["metric_name"]])) {
    scenario[["metric_name"]]
  } else if (identical(decision_direction, "gt")) {
    "non_rejection_rate"
  } else if (isTRUE(all.equal(as.numeric(scenario[["true_tau"]]), 0))) {
    "rejection_rate"
  } else {
    "power"
  }

  list(
    scenario_id = scenario[["id"]],
    n_simulations = n_target,
    n_attempts = attempts,
    metric_name = metric_name,
    metric_value = mean(metric_hits),
    metric_hit_count = as.integer(sum(metric_hits)),
    replication_trace = list(
      hit_indices = which(metric_hits),
      miss_indices = which(!metric_hits),
      hit_count = as.integer(sum(metric_hits)),
      miss_count = as.integer(sum(!metric_hits)),
      replications = replication_rows,
      treated_cluster_profile = treated_cluster_profile,
      treated_cluster_profile_summary = .summarize_treated_cluster_profile(
        treated_cluster_profile
      ),
      standard_t_counterfactual_summary =
        .summarize_story_local_wcb_standard_t_counterfactual(
          replication_rows = replication_rows,
          pvalue_threshold = pvalue_threshold,
          base_seed = scenario[["seed"]],
          required_rejection_count = as.integer(
            floor(as.numeric(scenario[["threshold"]]) * n_target) + 1L
          )
        ),
      non_near_threshold_miss_summary =
        .summarize_story_local_wcb_non_near_threshold_misses(
          replication_rows = replication_rows,
          pvalue_threshold = pvalue_threshold,
          base_seed = scenario[["seed"]],
          actual_rejection_count = as.integer(sum(metric_hits)),
          required_rejection_count = as.integer(
            floor(as.numeric(scenario[["threshold"]]) * n_target) + 1L
          )
        )
    ),
    metric_confidence_interval = .compute_exact_binomial_interval(
      successes = sum(metric_hits),
      trials = n_target
    ),
    requested_n_bootstrap = last_result$requested_n_bootstrap,
    actual_n_bootstrap = last_result$actual_n_bootstrap,
    weight_type = last_result$weight_type,
    impose_null = last_result$impose_null,
    full_enumeration = last_result$full_enumeration,
    ci_method = last_result$ci_method,
    pvalue_threshold = pvalue_threshold,
    decision_direction = decision_direction
  )
}


.validate_story_local_wcb_power_diagnostics_args <- function(
    threshold, true_tau, n_simulations, obs_per_cluster, seed, impose_null,
    replay_seeds = NULL,
    weight_type = "rademacher"
) {
  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      !is.finite(threshold) ||
      threshold < 0 ||
      threshold >= 1) {
    stop(
      "`threshold` must be a single finite number in [0, 1).",
      call. = FALSE
    )
  }

  if (!is.numeric(true_tau) || length(true_tau) != 1L || !is.finite(true_tau)) {
    stop("`true_tau` must be a single finite number.", call. = FALSE)
  }

  if (!is.numeric(n_simulations) ||
      length(n_simulations) != 1L ||
      !is.finite(n_simulations) ||
      n_simulations <= 0 ||
      n_simulations %% 1 != 0) {
    stop(
      "`n_simulations` must be a single positive integer.",
      call. = FALSE
    )
  }

  if (!is.numeric(obs_per_cluster) ||
      length(obs_per_cluster) != 1L ||
      !is.finite(obs_per_cluster) ||
      obs_per_cluster <= 0 ||
      obs_per_cluster %% 1 != 0) {
    stop(
      "`obs_per_cluster` must be a single positive integer.",
      call. = FALSE
    )
  }

  if (!is.numeric(seed) ||
      length(seed) != 1L ||
      !is.finite(seed) ||
      seed < 0 ||
      seed %% 1 != 0) {
    stop(
      "`seed` must be a single non-negative integer.",
      call. = FALSE
    )
  }

  if (!is.logical(impose_null) ||
      length(impose_null) != 1L ||
      is.na(impose_null)) {
    stop(
      "`impose_null` must be a single TRUE/FALSE value.",
      call. = FALSE
    )
  }

  supported_weight_types <- c("rademacher", "mammen", "webb")
  if (!is.character(weight_type) ||
      length(weight_type) != 1L ||
      is.na(weight_type) ||
      !weight_type %in% supported_weight_types) {
    stop(
      paste(
        "`weight_type` must be one of",
        "'rademacher', 'mammen', or 'webb'."
      ),
      call. = FALSE
    )
  }

  if (!is.null(replay_seeds)) {
    if (!is.numeric(replay_seeds) ||
        any(!is.finite(replay_seeds)) ||
        any(replay_seeds < 0) ||
        any(replay_seeds %% 1 != 0)) {
      stop(
        "`replay_seeds` must be NULL or a numeric vector of non-negative integers.",
        call. = FALSE
      )
    }
    if (any(as.integer(replay_seeds) < as.integer(seed))) {
      stop(
        paste(
          "`replay_seeds` must map to canonical sequential attempts",
          "at or after the scenario base seed."
        ),
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}


.replay_story_local_wcb_attempts <- function(
    scenario, shared_dgp = NULL, attempt_ids
) {
  attempt_ids <- as.integer(attempt_ids)
  if (length(attempt_ids) == 0L) {
    return(list())
  }
  if (any(attempt_ids < 1L)) {
    stop("Replay attempt ids must be positive integers.", call. = FALSE)
  }

  set.seed(as.integer(scenario[["seed"]]))

  alpha <- if (!is.null(scenario[["alpha"]])) {
    as.numeric(scenario[["alpha"]])
  } else {
    0.05
  }
  requested_n_bootstrap <- if (!is.null(scenario[["requested_n_bootstrap"]])) {
    as.integer(scenario[["requested_n_bootstrap"]])
  } else {
    as.integer(scenario[["n_bootstrap"]])
  }
  pvalue_threshold <- if (!is.null(scenario[["pvalue_threshold"]])) {
    as.numeric(scenario[["pvalue_threshold"]])
  } else {
    alpha
  }

  replay_rows <- vector("list", length(attempt_ids))
  max_attempt_id <- max(attempt_ids)
  attempts <- 0L
  valid_replications <- 0L

  while (attempts < max_attempt_id) {
    attempts <- attempts + 1L

    sim_data <- .simulate_clustering_monte_carlo_panel(scenario, shared_dgp)
    wild_result <- .run_wild_cluster_bootstrap_ci(
      data = sim_data,
      requested_n_bootstrap = requested_n_bootstrap,
      weight_type = scenario[["weight_type"]],
      alpha = alpha,
      impose_null = scenario[["impose_null"]]
    )

    if (!is.finite(wild_result$pvalue)) {
      next
    }

    valid_replications <- valid_replications + 1L
    match_idx <- which(attempt_ids == attempts)
    if (length(match_idx) == 0L) {
      next
    }

    replay_row <- list(
      attempt_id = as.integer(attempts),
      replication_id = as.integer(valid_replications),
      treated_cluster_count = as.integer(length(unique(sim_data$cluster[sim_data$D == 1]))),
      reject = as.logical(wild_result$pvalue < pvalue_threshold),
      pvalue = as.numeric(wild_result$pvalue),
      att = as.numeric(wild_result$att),
      original_se = as.numeric(wild_result$se),
      t_stat_original = as.numeric(wild_result$t_stat_original)
    )
    for (idx in match_idx) {
      replay_rows[[idx]] <- replay_row
    }
  }

  missing_rows <- which(vapply(replay_rows, is.null, logical(1L)))
  if (length(missing_rows) > 0L) {
    stop(
      sprintf(
        "Replay attempts %s did not yield finite wild-bootstrap p-values.",
        paste(attempt_ids[missing_rows], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  replay_rows
}


.build_story_local_wcb_targeted_replays <- function(
    scenario, shared_dgp, base_seed, replay_seeds, pvalue_threshold
) {
  if (is.null(replay_seeds)) {
    return(NULL)
  }

  replay_seeds <- as.integer(replay_seeds)
  replay_attempt_ids <- as.integer(replay_seeds - base_seed + 1L)

  replay_rows <- .replay_story_local_wcb_attempts(
    scenario = scenario,
    shared_dgp = shared_dgp,
    attempt_ids = replay_attempt_ids
  )

  lapply(seq_along(replay_rows), function(idx) {
    row <- replay_rows[[idx]]
    list(
      attempt_id = as.integer(row$attempt_id),
      replication_id = as.integer(row$replication_id),
      seed = as.integer(replay_seeds[[idx]]),
      replay_seed_offset = as.integer(replay_seeds[[idx]] - base_seed),
      treated_cluster_count = as.integer(row$treated_cluster_count),
      reject = isTRUE(row$reject),
      pvalue = as.numeric(row$pvalue),
      abs_gap_to_threshold = abs(as.numeric(row$pvalue) - pvalue_threshold),
      att = as.numeric(row$att),
      original_se = as.numeric(row$original_se),
      t_stat_original = as.numeric(row$t_stat_original)
    )
  })
}


.summarize_story_local_wcb_targeted_replays <- function(targeted_replays) {
  if (is.null(targeted_replays)) {
    return(NULL)
  }

  summarize_weakest_case <- function(reject_value) {
    matching_indices <- which(vapply(
      targeted_replays,
      function(row) identical(row$reject, reject_value),
      logical(1L)
    ))
    if (length(matching_indices) == 0L) {
      return(NULL)
    }

    candidate_gaps <- vapply(
      targeted_replays[matching_indices],
      `[[`,
      numeric(1L),
      "abs_gap_to_threshold"
    )
    weakest_row <- targeted_replays[[matching_indices[[which.min(candidate_gaps)]]]]

    list(
      seed = as.integer(weakest_row$seed),
      attempt_id = as.integer(weakest_row$attempt_id),
      replication_id = as.integer(weakest_row$replication_id),
      replay_seed_offset = as.integer(weakest_row$replay_seed_offset),
      treated_cluster_count = as.integer(weakest_row$treated_cluster_count),
      reject = isTRUE(weakest_row$reject),
      pvalue = as.numeric(weakest_row$pvalue),
      abs_gap_to_threshold = as.numeric(weakest_row$abs_gap_to_threshold)
    )
  }

  list(
    targeted_replay_count = as.integer(length(targeted_replays)),
    replay_seed_labels = vapply(
      targeted_replays,
      `[[`,
      integer(1L),
      "seed"
    ),
    replay_seed_offsets = vapply(
      targeted_replays,
      `[[`,
      integer(1L),
      "replay_seed_offset"
    ),
    replay_attempt_ids = vapply(
      targeted_replays,
      `[[`,
      integer(1L),
      "attempt_id"
    ),
    helper_rejection_pattern = vapply(
      targeted_replays,
      `[[`,
      logical(1L),
      "reject"
    ),
    treated_cluster_pattern = vapply(
      targeted_replays,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ),
    weakest_hit = summarize_weakest_case(TRUE),
    weakest_miss = summarize_weakest_case(FALSE)
  )
}


.compute_exact_binomial_interval <- function(
    successes, trials, conf_level = 0.95
) {
  interval <- stats::binom.test(
    x = as.integer(successes),
    n = as.integer(trials),
    conf.level = conf_level
  )$conf.int

  list(
    lower = unname(interval[[1L]]),
    upper = unname(interval[[2L]]),
    conf_level = as.numeric(conf_level)
  )
}


.summarize_story_local_wcb_bucket_repair_frontier <- function(
    bucket_ranking,
    actual_rejection_count,
    required_rejection_count,
    n_replications,
    non_near_threshold_miss_count
) {
  if (length(bucket_ranking) == 0L) {
    return(
      list(
        bucket_repair_frontier = list(),
        minimum_bucket_rank_to_clear_threshold = NULL,
        minimum_bucket_combo_to_clear_threshold = NULL,
        minimum_bucket_combo_share_of_non_near_threshold_misses = NULL
      )
    )
  }

  cumulative_miss_counts <- cumsum(vapply(
    bucket_ranking,
    `[[`,
    integer(1L),
    "miss_count"
  ))
  bucket_repair_frontier <- lapply(seq_along(bucket_ranking), function(idx) {
    bucket <- bucket_ranking[[idx]]
    cumulative_miss_count <- as.integer(cumulative_miss_counts[[idx]])
    counterfactual_rejection_count <- as.integer(
      actual_rejection_count + cumulative_miss_count
    )
    residual_shortfall <- as.integer(max(
      required_rejection_count - counterfactual_rejection_count,
      0L
    ))

    list(
      bucket_rank = as.integer(idx),
      treated_cluster_count = as.integer(bucket$treated_cluster_count),
      miss_count = as.integer(bucket$miss_count),
      cumulative_miss_count = cumulative_miss_count,
      counterfactual_rejection_count = counterfactual_rejection_count,
      counterfactual_power = counterfactual_rejection_count / n_replications,
      residual_shortfall = residual_shortfall,
      clears_threshold = identical(residual_shortfall, 0L),
      cumulative_share_of_non_near_threshold_misses = if (
        non_near_threshold_miss_count > 0L
      ) {
        cumulative_miss_count / non_near_threshold_miss_count
      } else {
        NaN
      }
    )
  })

  clear_indices <- which(vapply(
    bucket_repair_frontier,
    `[[`,
    logical(1L),
    "clears_threshold"
  ))
  if (length(clear_indices) == 0L) {
    return(
      list(
        bucket_repair_frontier = bucket_repair_frontier,
        minimum_bucket_rank_to_clear_threshold = NULL,
        minimum_bucket_combo_to_clear_threshold = NULL,
        minimum_bucket_combo_share_of_non_near_threshold_misses = NULL
      )
    )
  }

  first_clear_idx <- clear_indices[[1L]]

  list(
    bucket_repair_frontier = bucket_repair_frontier,
    minimum_bucket_rank_to_clear_threshold = as.integer(first_clear_idx),
    minimum_bucket_combo_to_clear_threshold = vapply(
      bucket_ranking[seq_len(first_clear_idx)],
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ),
    minimum_bucket_combo_share_of_non_near_threshold_misses =
      bucket_repair_frontier[[first_clear_idx]]$cumulative_share_of_non_near_threshold_misses
  )
}


.summarize_story_local_wcb_non_near_threshold_case_repair_frontier <- function(
    rows,
    base_seed,
    actual_rejection_count,
    required_rejection_count,
    total_replication_count,
    pvalue_threshold
) {
  non_near_threshold_miss_count <- as.integer(length(rows))
  if (non_near_threshold_miss_count == 0L) {
    return(
      list(
        case_repair_frontier = list(),
        minimum_case_rank_to_clear_threshold = NULL,
        minimum_case_replay_seed_combo_to_clear_threshold = NULL,
        minimum_case_combo_share_of_non_near_threshold_misses = NULL
      )
    )
  }

  ordered_rows <- rows[order(
    vapply(rows, `[[`, numeric(1L), "pvalue"),
    vapply(rows, `[[`, integer(1L), "attempt_id")
  )]
  frontier <- lapply(seq_along(ordered_rows), function(idx) {
    row <- ordered_rows[[idx]]
    replay_seed <- as.integer(base_seed + row$attempt_id - 1L)
    cumulative_case_count <- as.integer(idx)
    counterfactual_rejection_count <- as.integer(
      actual_rejection_count + cumulative_case_count
    )
    residual_shortfall <- as.integer(max(
      required_rejection_count - counterfactual_rejection_count,
      0L
    ))

    list(
      case_rank = as.integer(idx),
      attempt_id = as.integer(row$attempt_id),
      replication_id = as.integer(row$replication_id),
      replay_seed = replay_seed,
      replay_seed_offset = as.integer(row$attempt_id - 1L),
      treated_cluster_count = as.integer(row$treated_cluster_count),
      pvalue = as.numeric(row$pvalue),
      pvalue_gap_above_threshold = as.numeric(row$pvalue - pvalue_threshold),
      cumulative_case_count = cumulative_case_count,
      counterfactual_rejection_count = counterfactual_rejection_count,
      counterfactual_power = counterfactual_rejection_count / total_replication_count,
      residual_shortfall = residual_shortfall,
      clears_threshold = identical(residual_shortfall, 0L),
      cumulative_share_of_non_near_threshold_misses =
        cumulative_case_count / non_near_threshold_miss_count
    )
  })

  clear_indices <- which(vapply(
    frontier,
    `[[`,
    logical(1L),
    "clears_threshold"
  ))
  if (length(clear_indices) == 0L) {
    return(
      list(
        case_repair_frontier = frontier,
        minimum_case_rank_to_clear_threshold = NULL,
        minimum_case_replay_seed_combo_to_clear_threshold = NULL,
        minimum_case_combo_share_of_non_near_threshold_misses = NULL
      )
    )
  }

  first_clear_idx <- clear_indices[[1L]]

  list(
    case_repair_frontier = frontier,
    minimum_case_rank_to_clear_threshold = as.integer(first_clear_idx),
    minimum_case_replay_seed_combo_to_clear_threshold = vapply(
      frontier[seq_len(first_clear_idx)],
      `[[`,
      integer(1L),
      "replay_seed"
    ),
    minimum_case_combo_share_of_non_near_threshold_misses =
      frontier[[first_clear_idx]]$cumulative_share_of_non_near_threshold_misses
  )
}


.classify_story_local_wcb_representative_far_miss_mechanism <- function(
    representative_large_miss
) {
  low_signal <- is.finite(
    representative_large_miss$abs_att_to_bucket_avg_ratio
  ) && representative_large_miss$abs_att_to_bucket_avg_ratio < 1
  elevated_se <- is.finite(
    representative_large_miss$original_se_to_bucket_avg_ratio
  ) && representative_large_miss$original_se_to_bucket_avg_ratio > 1

  if (low_signal && elevated_se) {
    return("low-signal-plus-elevated-se")
  }

  if (low_signal) {
    return("low-signal-dominant")
  }

  if (elevated_se) {
    return("elevated-se-dominant")
  }

  "mixed"
}


.summarize_story_local_wcb_representative_far_miss_mechanism <- function(
    bucket_ranking
) {
  if (length(bucket_ranking) == 0L) {
    return(
      list(
        top_bucket_order = integer(0L),
        representative_replay_seeds = integer(0L),
        low_signal_plus_elevated_se_buckets = integer(0L),
        low_signal_dominant_buckets = integer(0L),
        mechanism_split = NULL,
        all_representative_abs_att_ratios_below_one = NA,
        all_buckets_ranked_by_non_near_threshold_mass = TRUE,
        bucket_mechanisms = list()
      )
    )
  }

  top_bucket_limit <- min(length(bucket_ranking), 3L)
  top_buckets <- bucket_ranking[seq_len(top_bucket_limit)]

  bucket_mechanisms <- lapply(top_buckets, function(bucket) {
    representative_large_miss <- bucket$representative_large_miss
    mechanism_label <- .classify_story_local_wcb_representative_far_miss_mechanism(
      representative_large_miss = representative_large_miss
    )

    list(
      treated_cluster_count = as.integer(bucket$treated_cluster_count),
      representative_replay_seed = as.integer(representative_large_miss$replay_seed),
      mechanism_label = mechanism_label,
      representative_abs_att_to_bucket_avg_ratio = as.numeric(
        representative_large_miss$abs_att_to_bucket_avg_ratio
      ),
      representative_original_se_to_bucket_avg_ratio = as.numeric(
        representative_large_miss$original_se_to_bucket_avg_ratio
      )
    )
  })

  top_bucket_order <- as.integer(vapply(
    bucket_mechanisms,
    `[[`,
    integer(1L),
    "treated_cluster_count"
  ))
  representative_replay_seeds <- as.integer(vapply(
    bucket_mechanisms,
    `[[`,
    integer(1L),
    "representative_replay_seed"
  ))
  mechanism_labels <- vapply(
    bucket_mechanisms,
    `[[`,
    character(1L),
    "mechanism_label"
  )
  low_signal_plus_elevated_se_buckets <- top_bucket_order[
    mechanism_labels == "low-signal-plus-elevated-se"
  ]
  low_signal_dominant_buckets <- top_bucket_order[
    mechanism_labels == "low-signal-dominant"
  ]
  split_parts <- character(0L)
  if (length(low_signal_plus_elevated_se_buckets) > 0L) {
    split_parts <- c(
      split_parts,
      sprintf(
        "buckets-%s-low-signal-plus-elevated-se",
        paste(low_signal_plus_elevated_se_buckets, collapse = "-")
      )
    )
  }
  if (length(low_signal_dominant_buckets) > 0L) {
    split_parts <- c(
      split_parts,
      sprintf(
        "bucket-%s-low-signal-dominant",
        paste(low_signal_dominant_buckets, collapse = "-")
      )
    )
  }
  mechanism_split <- if (length(split_parts) > 0L) {
    paste(split_parts, collapse = "_")
  } else {
    NULL
  }

  list(
    top_bucket_order = top_bucket_order,
    representative_replay_seeds = representative_replay_seeds,
    low_signal_plus_elevated_se_buckets = as.integer(
      low_signal_plus_elevated_se_buckets
    ),
    low_signal_dominant_buckets = as.integer(low_signal_dominant_buckets),
    mechanism_split = mechanism_split,
    all_representative_abs_att_ratios_below_one = all(vapply(
      bucket_mechanisms,
      function(bucket) {
        is.finite(bucket$representative_abs_att_to_bucket_avg_ratio) &&
          bucket$representative_abs_att_to_bucket_avg_ratio < 1
      },
      logical(1L)
    )),
    all_buckets_ranked_by_non_near_threshold_mass = identical(
      top_bucket_order,
      as.integer(vapply(
        top_buckets,
        `[[`,
        integer(1L),
        "treated_cluster_count"
      ))
    ),
    bucket_mechanisms = bucket_mechanisms
  )
}


.summarize_story_local_wcb_non_near_threshold_misses <- function(
    replication_rows,
    pvalue_threshold,
    base_seed,
    actual_rejection_count,
    required_rejection_count
) {
  base_seed <- as.integer(base_seed)
  near_threshold_window <- list(
    lower = as.numeric(round(pvalue_threshold - 0.01, 2)),
    upper = as.numeric(round(pvalue_threshold + 0.01, 2))
  )

  if (length(replication_rows) == 0L) {
    return(
      list(
        near_threshold_window = near_threshold_window,
        total_misses = 0L,
        near_threshold_misses = 0L,
        non_near_threshold_misses = 0L,
        non_near_threshold_share_of_misses = NaN,
        bucket_ranking = list(),
        dominant_bucket = NULL,
        bucket_repair_frontier = list(),
        minimum_bucket_rank_to_clear_threshold = NULL,
        minimum_bucket_combo_to_clear_threshold = NULL,
        minimum_bucket_combo_share_of_non_near_threshold_misses = NULL,
        case_repair_frontier = list(),
        minimum_case_rank_to_clear_threshold = NULL,
        minimum_case_replay_seed_combo_to_clear_threshold = NULL,
        minimum_case_combo_share_of_non_near_threshold_misses = NULL,
        representative_far_miss_mechanism_summary =
          .summarize_story_local_wcb_representative_far_miss_mechanism(
            bucket_ranking = list()
          )
      )
    )
  }

  miss_rows <- replication_rows[!vapply(
    replication_rows,
    `[[`,
    logical(1L),
    "reject"
  )]
  near_threshold_miss_rows <- miss_rows[vapply(
    miss_rows,
    function(row) {
      row$pvalue >= near_threshold_window$lower &&
        row$pvalue <= near_threshold_window$upper
    },
    logical(1L)
  )]
  non_near_threshold_miss_rows <- miss_rows[vapply(
    miss_rows,
    function(row) row$pvalue > near_threshold_window$upper,
    logical(1L)
  )]

  non_near_threshold_count <- length(non_near_threshold_miss_rows)
  bucket_ranking <- list()
  dominant_bucket <- NULL

  if (non_near_threshold_count > 0L) {
    bucket_names <- vapply(
      non_near_threshold_miss_rows,
      function(row) as.character(row$treated_cluster_count),
      character(1L)
    )
    bucket_rows <- split(non_near_threshold_miss_rows, bucket_names)
    ordered_bucket_names <- names(bucket_rows)[order(
      -vapply(bucket_rows, length, integer(1L)),
      as.integer(names(bucket_rows))
    )]

    bucket_ranking <- lapply(ordered_bucket_names, function(bucket_name) {
      rows <- bucket_rows[[bucket_name]]
      bucket_pvalues <- vapply(rows, `[[`, numeric(1L), "pvalue")
      bucket_abs_att <- abs(vapply(
        rows,
        `[[`,
        numeric(1L),
        "att"
      ))
      bucket_original_se <- vapply(
        rows,
        `[[`,
        numeric(1L),
        "original_se"
      )
      bucket_abs_t_stats <- abs(vapply(
        rows,
        `[[`,
        numeric(1L),
        "t_stat_original"
      ))
      bucket_attempt_ids <- vapply(rows, `[[`, integer(1L), "attempt_id")
      bucket_replay_seeds <- base_seed + bucket_attempt_ids - 1L
      representative_idx <- order(-bucket_pvalues, bucket_replay_seeds)[[1L]]

      list(
        treated_cluster_count = as.integer(bucket_name),
        miss_count = as.integer(length(rows)),
        share_of_non_near_threshold_misses = length(rows) / non_near_threshold_count,
        avg_pvalue = mean(bucket_pvalues),
        median_pvalue = stats::median(bucket_pvalues),
        avg_abs_att = mean(bucket_abs_att),
        median_abs_att = stats::median(bucket_abs_att),
        avg_original_se = mean(bucket_original_se),
        median_original_se = stats::median(bucket_original_se),
        max_pvalue = max(bucket_pvalues),
        avg_abs_t_stat = mean(bucket_abs_t_stats),
        median_abs_t_stat = stats::median(bucket_abs_t_stats),
        representative_large_miss = list(
          attempt_id = as.integer(bucket_attempt_ids[[representative_idx]]),
          replay_seed = as.integer(bucket_replay_seeds[[representative_idx]]),
          replay_seed_offset = as.integer(
            bucket_replay_seeds[[representative_idx]] - base_seed
          ),
          pvalue = as.numeric(bucket_pvalues[[representative_idx]]),
          pvalue_gap = as.numeric(
            bucket_pvalues[[representative_idx]] - pvalue_threshold
          ),
          abs_att = as.numeric(bucket_abs_att[[representative_idx]]),
          original_se = as.numeric(bucket_original_se[[representative_idx]]),
          abs_att_to_bucket_avg_ratio = as.numeric(
            bucket_abs_att[[representative_idx]] / mean(bucket_abs_att)
          ),
          original_se_to_bucket_avg_ratio = as.numeric(
            bucket_original_se[[representative_idx]] / mean(bucket_original_se)
          ),
          abs_t_stat_to_bucket_avg_ratio = as.numeric(
            bucket_abs_t_stats[[representative_idx]] / mean(bucket_abs_t_stats)
          ),
          abs_t_stat_original = as.numeric(
            bucket_abs_t_stats[[representative_idx]]
          )
        )
      )
    })
    dominant_bucket <- bucket_ranking[[1L]]
  }

  total_misses <- as.integer(length(miss_rows))
  near_threshold_misses <- as.integer(length(near_threshold_miss_rows))
  non_near_threshold_misses <- as.integer(non_near_threshold_count)
  bucket_repair_summary <- .summarize_story_local_wcb_bucket_repair_frontier(
    bucket_ranking = bucket_ranking,
    actual_rejection_count = as.integer(actual_rejection_count),
    required_rejection_count = as.integer(required_rejection_count),
    n_replications = as.integer(length(replication_rows)),
    non_near_threshold_miss_count = non_near_threshold_misses
  )
  case_repair_summary <-
    .summarize_story_local_wcb_non_near_threshold_case_repair_frontier(
      rows = non_near_threshold_miss_rows,
      base_seed = base_seed,
      actual_rejection_count = as.integer(actual_rejection_count),
      required_rejection_count = as.integer(required_rejection_count),
      total_replication_count = as.integer(length(replication_rows)),
      pvalue_threshold = pvalue_threshold
    )
  list(
    near_threshold_window = near_threshold_window,
    total_misses = total_misses,
    near_threshold_misses = near_threshold_misses,
    non_near_threshold_misses = non_near_threshold_misses,
    non_near_threshold_share_of_misses = if (total_misses > 0L) {
      non_near_threshold_misses / total_misses
    } else {
      NaN
    },
    bucket_ranking = bucket_ranking,
    dominant_bucket = dominant_bucket,
    bucket_repair_frontier = bucket_repair_summary$bucket_repair_frontier,
    minimum_bucket_rank_to_clear_threshold =
      bucket_repair_summary$minimum_bucket_rank_to_clear_threshold,
    minimum_bucket_combo_to_clear_threshold =
      bucket_repair_summary$minimum_bucket_combo_to_clear_threshold,
    minimum_bucket_combo_share_of_non_near_threshold_misses =
      bucket_repair_summary$minimum_bucket_combo_share_of_non_near_threshold_misses,
    case_repair_frontier = case_repair_summary$case_repair_frontier,
    minimum_case_rank_to_clear_threshold =
      case_repair_summary$minimum_case_rank_to_clear_threshold,
    minimum_case_replay_seed_combo_to_clear_threshold =
      case_repair_summary$minimum_case_replay_seed_combo_to_clear_threshold,
    minimum_case_combo_share_of_non_near_threshold_misses =
      case_repair_summary$minimum_case_combo_share_of_non_near_threshold_misses,
    representative_far_miss_mechanism_summary =
      .summarize_story_local_wcb_representative_far_miss_mechanism(
        bucket_ranking = bucket_ranking
      )
  )
}


.summarize_story_local_wcb_standard_t_counterfactual <- function(
    replication_rows,
    pvalue_threshold,
    base_seed,
    required_rejection_count
) {
  total_replication_count <- as.integer(length(replication_rows))
  required_rejection_count <- as.integer(required_rejection_count)
  near_threshold_window <- list(
    lower = as.numeric(round(pvalue_threshold - 0.01, 2)),
    upper = as.numeric(round(pvalue_threshold + 0.01, 2))
  )

  build_counterfactual_cases <- function(rows) {
    if (length(rows) == 0L) {
      return(list())
    }

    lapply(rows, function(row) {
      list(
        attempt_id = as.integer(row$attempt_id),
        replication_id = as.integer(row$replication_id),
        replay_seed = as.integer(base_seed + row$attempt_id - 1L),
        replay_seed_offset = as.integer(row$attempt_id - 1L),
        treated_cluster_count = as.integer(row$treated_cluster_count),
        wcb_pvalue = as.numeric(row$pvalue),
        standard_t_pvalue = as.numeric(row$standard_t_pvalue)
      )
    })
  }

  summarize_buckets <- function(rows) {
    if (length(rows) == 0L) {
      return(list())
    }

    bucket_counts <- table(vapply(
      rows,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ))
    ordered_names <- names(bucket_counts)[order(
      -as.integer(bucket_counts),
      as.integer(names(bucket_counts))
    )]

    lapply(ordered_names, function(bucket_name) {
      list(
        treated_cluster_count = as.integer(bucket_name),
        count = as.integer(bucket_counts[[bucket_name]])
      )
    })
  }

  extract_replay_seeds <- function(rows) {
    if (length(rows) == 0L) {
      return(integer(0L))
    }

    as.integer(base_seed + vapply(
      rows,
      `[[`,
      integer(1L),
      "attempt_id"
    ) - 1L)
  }

  summarize_closest_case <- function(rows, boundary) {
    if (length(rows) == 0L) {
      return(list(seed = NULL, gap = NULL))
    }

    gaps <- vapply(
      rows,
      function(row) as.numeric(row$pvalue - boundary),
      numeric(1L)
    )
    closest_index <- which.min(gaps)
    closest_row <- rows[[closest_index]]

    list(
      seed = as.integer(base_seed + closest_row$attempt_id - 1L),
      gap = gaps[[closest_index]]
    )
  }

  summarize_case_repair_frontier <- function(
      rows,
      actual_rejection_count,
      required_rejection_count
  ) {
    standard_only_rejection_count <- as.integer(length(rows))
    if (standard_only_rejection_count == 0L) {
      return(
        list(
          standard_only_repair_frontier = list(),
          minimum_case_rank_to_clear_threshold = NULL,
          minimum_case_replay_seed_combo_to_clear_threshold = NULL,
          minimum_case_combo_share_of_standard_only_cases = NULL
        )
      )
    }

    ordered_rows <- rows[order(
      vapply(rows, `[[`, numeric(1L), "pvalue"),
      vapply(rows, `[[`, integer(1L), "attempt_id")
    )]
    frontier <- lapply(seq_along(ordered_rows), function(idx) {
      row <- ordered_rows[[idx]]
      replay_seed <- as.integer(base_seed + row$attempt_id - 1L)
      cumulative_rescue_count <- as.integer(idx)
      counterfactual_rejection_count <- as.integer(
        actual_rejection_count + cumulative_rescue_count
      )
      residual_shortfall <- as.integer(max(
        required_rejection_count - counterfactual_rejection_count,
        0L
      ))

      list(
        case_rank = as.integer(idx),
        attempt_id = as.integer(row$attempt_id),
        replication_id = as.integer(row$replication_id),
        replay_seed = replay_seed,
        replay_seed_offset = as.integer(row$attempt_id - 1L),
        treated_cluster_count = as.integer(row$treated_cluster_count),
        wcb_pvalue = as.numeric(row$pvalue),
        standard_t_pvalue = as.numeric(row$standard_t_pvalue),
        wcb_pvalue_gap = as.numeric(row$pvalue - pvalue_threshold),
        paired_pvalue_gap = as.numeric(row$pvalue - row$standard_t_pvalue),
        cumulative_rescue_count = cumulative_rescue_count,
        counterfactual_rejection_count = counterfactual_rejection_count,
        counterfactual_power = counterfactual_rejection_count / total_replication_count,
        residual_shortfall = residual_shortfall,
        clears_threshold = identical(residual_shortfall, 0L),
        cumulative_share_of_standard_only_cases =
          cumulative_rescue_count / standard_only_rejection_count
      )
    })

    clear_indices <- which(vapply(
      frontier,
      `[[`,
      logical(1L),
      "clears_threshold"
    ))
    if (length(clear_indices) == 0L) {
      return(
        list(
          standard_only_repair_frontier = frontier,
          minimum_case_rank_to_clear_threshold = NULL,
          minimum_case_replay_seed_combo_to_clear_threshold = NULL,
          minimum_case_combo_share_of_standard_only_cases = NULL
        )
      )
    }

    first_clear_idx <- clear_indices[[1L]]

    list(
      standard_only_repair_frontier = frontier,
      minimum_case_rank_to_clear_threshold = as.integer(first_clear_idx),
      minimum_case_replay_seed_combo_to_clear_threshold = vapply(
        frontier[seq_len(first_clear_idx)],
        `[[`,
        integer(1L),
        "replay_seed"
      ),
      minimum_case_combo_share_of_standard_only_cases =
        frontier[[first_clear_idx]]$cumulative_share_of_standard_only_cases
    )
  }

  paired_rows <- replication_rows[vapply(
    replication_rows,
    function(row) {
      is.logical(row$standard_t_reject) &&
        length(row$standard_t_reject) == 1L &&
        !is.na(row$standard_t_reject) &&
        is.finite(row$standard_t_pvalue)
    },
    logical(1L)
  )]
  if (length(paired_rows) == 0L) {
    return(
      list(
        paired_replication_count = 0L,
        standard_t_rejection_count = 0L,
        wild_bootstrap_rejection_count = 0L,
        shared_rejection_count = 0L,
        standard_only_rejection_count = 0L,
        wcb_only_rejection_count = 0L,
        paired_rejection_gap = 0L,
        counterfactual_power = if (total_replication_count > 0L) 0 else NaN,
        required_rejection_count_to_clear_threshold = required_rejection_count,
        required_power_to_clear_threshold = if (total_replication_count > 0L) {
          required_rejection_count / total_replication_count
        } else {
          NaN
        },
        residual_rejection_shortfall = required_rejection_count,
        counterfactual_clears_threshold = identical(required_rejection_count, 0L),
        near_threshold_window = near_threshold_window,
        standard_only_near_threshold_count = 0L,
        standard_only_non_near_threshold_count = 0L,
        standard_only_near_threshold_share = NaN,
        standard_only_non_near_threshold_share = NaN,
        standard_only_near_threshold_replication_share = NaN,
        standard_only_non_near_threshold_replication_share = NaN,
        standard_only_near_threshold_replay_seeds = integer(0L),
        standard_only_non_near_threshold_replay_seeds = integer(0L),
        closest_standard_only_seed_to_threshold = NULL,
        closest_standard_only_gap_to_threshold = NULL,
        closest_non_near_standard_only_seed = NULL,
        closest_non_near_standard_only_gap_to_window_upper = NULL,
        standard_only_near_threshold_bucket_ranking = list(),
        standard_only_non_near_threshold_bucket_ranking = list(),
        standard_only_replay_seeds = integer(0L),
        standard_only_bucket_ranking = list(),
        standard_only_cases = list(),
        wcb_only_cases = list(),
        standard_only_repair_frontier = list(),
        minimum_case_rank_to_clear_threshold = NULL,
        minimum_case_replay_seed_combo_to_clear_threshold = NULL,
        minimum_case_combo_share_of_standard_only_cases = NULL
      )
    )
  }

  standard_t_reject <- vapply(
    paired_rows,
    `[[`,
    logical(1L),
    "standard_t_reject"
  )
  wild_bootstrap_reject <- vapply(
    paired_rows,
    `[[`,
    logical(1L),
    "reject"
  )
  standard_only_rows <- paired_rows[standard_t_reject & !wild_bootstrap_reject]
  wcb_only_rows <- paired_rows[!standard_t_reject & wild_bootstrap_reject]
  standard_only_near_threshold_rows <- standard_only_rows[vapply(
    standard_only_rows,
    function(row) {
      row$pvalue >= near_threshold_window$lower &&
        row$pvalue <= near_threshold_window$upper
    },
    logical(1L)
  )]
  standard_only_non_near_threshold_rows <- standard_only_rows[vapply(
    standard_only_rows,
    function(row) row$pvalue > near_threshold_window$upper,
    logical(1L)
  )]
  standard_only_cases <- build_counterfactual_cases(standard_only_rows)
  wcb_only_cases <- build_counterfactual_cases(wcb_only_rows)
  standard_only_rejection_count <- as.integer(length(standard_only_rows))
  closest_standard_only_row <- if (standard_only_rejection_count > 0L) {
    standard_only_rows[[order(
      vapply(
        standard_only_rows,
        function(row) row$pvalue - pvalue_threshold,
        numeric(1L)
      ),
      vapply(standard_only_rows, `[[`, integer(1L), "attempt_id")
    )[[1L]]]]
  } else {
    NULL
  }
  closest_non_near_standard_only_row <- if (
    length(standard_only_non_near_threshold_rows) > 0L
  ) {
    standard_only_non_near_threshold_rows[[order(
      vapply(
        standard_only_non_near_threshold_rows,
        function(row) row$pvalue - near_threshold_window$upper,
        numeric(1L)
      ),
      vapply(
        standard_only_non_near_threshold_rows,
        `[[`,
        integer(1L),
        "attempt_id"
      )
    )[[1L]]]]
  } else {
    NULL
  }
  standard_t_rejection_count <- as.integer(sum(standard_t_reject))
  closest_near_threshold_case <- summarize_closest_case(
    rows = standard_only_near_threshold_rows,
    boundary = pvalue_threshold
  )
  closest_non_near_threshold_case <- summarize_closest_case(
    rows = standard_only_non_near_threshold_rows,
    boundary = near_threshold_window$upper
  )
  case_repair_summary <- summarize_case_repair_frontier(
    rows = standard_only_rows,
    actual_rejection_count = as.integer(sum(wild_bootstrap_reject)),
    required_rejection_count = required_rejection_count
  )

  list(
    paired_replication_count = as.integer(length(paired_rows)),
    standard_t_rejection_count = standard_t_rejection_count,
    wild_bootstrap_rejection_count = as.integer(sum(wild_bootstrap_reject)),
    shared_rejection_count = as.integer(sum(
      standard_t_reject & wild_bootstrap_reject
    )),
    standard_only_rejection_count = as.integer(length(standard_only_rows)),
    wcb_only_rejection_count = as.integer(length(wcb_only_rows)),
    paired_rejection_gap = as.integer(sum(standard_t_reject) - sum(wild_bootstrap_reject)),
    counterfactual_power = if (total_replication_count > 0L) {
      standard_t_rejection_count / total_replication_count
    } else {
      NaN
    },
    required_rejection_count_to_clear_threshold = required_rejection_count,
    required_power_to_clear_threshold = if (total_replication_count > 0L) {
      required_rejection_count / total_replication_count
    } else {
      NaN
    },
    residual_rejection_shortfall = as.integer(max(
      required_rejection_count - standard_t_rejection_count,
      0L
    )),
    counterfactual_clears_threshold = standard_t_rejection_count >=
      required_rejection_count,
    near_threshold_window = near_threshold_window,
    standard_only_near_threshold_count = as.integer(length(
      standard_only_near_threshold_rows
    )),
    standard_only_non_near_threshold_count = as.integer(length(
      standard_only_non_near_threshold_rows
    )),
    standard_only_near_threshold_share = if (length(standard_only_rows) > 0L) {
      length(standard_only_near_threshold_rows) / length(standard_only_rows)
    } else {
      NaN
    },
    standard_only_non_near_threshold_share = if (length(standard_only_rows) > 0L) {
      length(standard_only_non_near_threshold_rows) / length(standard_only_rows)
    } else {
      NaN
    },
    standard_only_near_threshold_replication_share = if (total_replication_count > 0L) {
      length(standard_only_near_threshold_rows) / total_replication_count
    } else {
      NaN
    },
    standard_only_non_near_threshold_replication_share = if (
      total_replication_count > 0L
    ) {
      length(standard_only_non_near_threshold_rows) / total_replication_count
    } else {
      NaN
    },
    standard_only_near_threshold_replay_seeds = extract_replay_seeds(
      standard_only_near_threshold_rows
    ),
    standard_only_non_near_threshold_replay_seeds = extract_replay_seeds(
      standard_only_non_near_threshold_rows
    ),
    closest_standard_only_seed_to_threshold = closest_near_threshold_case$seed,
    closest_standard_only_gap_to_threshold = closest_near_threshold_case$gap,
    closest_non_near_standard_only_seed =
      closest_non_near_threshold_case$seed,
    closest_non_near_standard_only_gap_to_window_upper =
      closest_non_near_threshold_case$gap,
    standard_only_near_threshold_bucket_ranking = summarize_buckets(
      standard_only_near_threshold_rows
    ),
    standard_only_non_near_threshold_bucket_ranking = summarize_buckets(
      standard_only_non_near_threshold_rows
    ),
    standard_only_replay_seeds = as.integer(vapply(
      standard_only_cases,
      `[[`,
      integer(1L),
      "replay_seed"
    )),
    standard_only_bucket_ranking = summarize_buckets(standard_only_rows),
    standard_only_cases = standard_only_cases,
    wcb_only_cases = wcb_only_cases,
    standard_only_repair_frontier =
      case_repair_summary$standard_only_repair_frontier,
    minimum_case_rank_to_clear_threshold =
      case_repair_summary$minimum_case_rank_to_clear_threshold,
    minimum_case_replay_seed_combo_to_clear_threshold =
      case_repair_summary$minimum_case_replay_seed_combo_to_clear_threshold,
    minimum_case_combo_share_of_standard_only_cases =
      case_repair_summary$minimum_case_combo_share_of_standard_only_cases
  )
}


.summarize_treated_cluster_profile <- function(treated_cluster_profile) {
  summarize_bucket_group <- function(buckets) {
    if (length(buckets) == 0L) {
      return(
        list(
          treated_cluster_counts = integer(0L),
          total_replications = 0L,
          total_rejections = 0L,
          share_of_rejections = NaN,
          share_of_replications = NaN
        )
      )
    }

    total_replications <- sum(vapply(
      buckets,
      `[[`,
      integer(1L),
      "n_replications"
    ))
    total_rejections <- sum(vapply(
      buckets,
      `[[`,
      integer(1L),
      "rejection_count"
    ))

    list(
      treated_cluster_counts = vapply(
        buckets,
        `[[`,
        integer(1L),
        "treated_cluster_count"
      ),
      total_replications = as.integer(total_replications),
      total_rejections = as.integer(total_rejections),
      share_of_rejections = NaN,
      share_of_replications = NaN
    )
  }

  if (length(treated_cluster_profile) == 0L) {
    zero_summary <- summarize_bucket_group(list())

    return(
      list(
        modal_bucket = NULL,
        modal_bucket_signal_ratio_contrast = NULL,
        high_rejection_buckets = zero_summary,
        zero_rejection_buckets = list(
          treated_cluster_counts = zero_summary$treated_cluster_counts,
          total_replications = zero_summary$total_replications
        ),
        near_threshold_buckets = list(
          treated_cluster_counts = zero_summary$treated_cluster_counts,
          total_near_threshold_cases = 0L
        )
      )
    )
  }

  bucket_replications <- vapply(
    treated_cluster_profile,
    `[[`,
    integer(1L),
    "n_replications"
  )
  bucket_rates <- vapply(
    treated_cluster_profile,
    `[[`,
    numeric(1L),
    "rejection_rate"
  )
  bucket_near_threshold <- vapply(
    treated_cluster_profile,
    `[[`,
    integer(1L),
    "near_threshold_count"
  )
  total_replications <- sum(bucket_replications)
  total_rejections <- sum(vapply(
    treated_cluster_profile,
    `[[`,
    integer(1L),
    "rejection_count"
  ))

  high_rejection_summary <- summarize_bucket_group(
    treated_cluster_profile[bucket_rates >= 0.5]
  )
  high_rejection_summary$share_of_rejections <- if (total_rejections > 0L) {
    high_rejection_summary$total_rejections / total_rejections
  } else {
    NaN
  }
  high_rejection_summary$share_of_replications <- if (total_replications > 0L) {
    high_rejection_summary$total_replications / total_replications
  } else {
    NaN
  }

  high_rejection_buckets <- treated_cluster_profile[bucket_rates >= 0.5]
  modal_bucket <- treated_cluster_profile[[which.max(bucket_replications)]]
  comparison_candidates <- high_rejection_buckets[vapply(
    high_rejection_buckets,
    function(bucket) bucket$treated_cluster_count != modal_bucket$treated_cluster_count,
    logical(1L)
  )]
  comparison_bucket_is_high_rejection <- TRUE
  if (length(comparison_candidates) == 0L) {
    comparison_candidates <- treated_cluster_profile[vapply(
      treated_cluster_profile,
      function(bucket) bucket$treated_cluster_count != modal_bucket$treated_cluster_count,
      logical(1L)
    )]
    comparison_bucket_is_high_rejection <- FALSE
  }
  modal_bucket_signal_ratio_contrast <- if (length(comparison_candidates) == 0L) {
    NULL
  } else {
    comparison_bucket <- comparison_candidates[[which.max(vapply(
      comparison_candidates,
      `[[`,
      numeric(1L),
      "avg_abs_att_to_se_ratio"
    ))]]
    avg_ratio_gap <- as.numeric(
      modal_bucket$avg_abs_att_to_se_ratio - comparison_bucket$avg_abs_att_to_se_ratio
    )
    median_ratio_gap <- as.numeric(
      modal_bucket$median_abs_att_to_se_ratio -
        comparison_bucket$median_abs_att_to_se_ratio
    )

    list(
      modal_bucket_treated_cluster_count = as.integer(modal_bucket$treated_cluster_count),
      comparison_bucket_treated_cluster_count = as.integer(
        comparison_bucket$treated_cluster_count
      ),
      comparison_bucket_is_high_rejection = comparison_bucket_is_high_rejection,
      modal_avg_abs_att_to_se_ratio = as.numeric(modal_bucket$avg_abs_att_to_se_ratio),
      comparison_avg_abs_att_to_se_ratio = as.numeric(
        comparison_bucket$avg_abs_att_to_se_ratio
      ),
      modal_median_abs_att_to_se_ratio = as.numeric(
        modal_bucket$median_abs_att_to_se_ratio
      ),
      comparison_median_abs_att_to_se_ratio = as.numeric(
        comparison_bucket$median_abs_att_to_se_ratio
      ),
      avg_ratio_gap_modal_minus_comparison = avg_ratio_gap,
      median_ratio_gap_modal_minus_comparison = median_ratio_gap,
      modal_bucket_weaker_on_avg_ratio = isTRUE(avg_ratio_gap < 0),
      modal_bucket_weaker_on_median_ratio = isTRUE(median_ratio_gap < 0)
    )
  }

  zero_rejection_summary <- summarize_bucket_group(
    treated_cluster_profile[vapply(
      treated_cluster_profile,
      function(bucket) bucket$rejection_count == 0L,
      logical(1L)
    )]
  )
  near_threshold_buckets <- treated_cluster_profile[bucket_near_threshold > 0L]
  near_threshold_summary <- summarize_bucket_group(near_threshold_buckets)

  list(
    modal_bucket = modal_bucket,
    modal_bucket_signal_ratio_contrast = modal_bucket_signal_ratio_contrast,
    high_rejection_buckets = high_rejection_summary,
    zero_rejection_buckets = list(
      treated_cluster_counts = zero_rejection_summary$treated_cluster_counts,
      total_replications = zero_rejection_summary$total_replications
    ),
    near_threshold_buckets = list(
      treated_cluster_counts = near_threshold_summary$treated_cluster_counts,
      total_near_threshold_cases = as.integer(sum(vapply(
        near_threshold_buckets,
        `[[`,
        integer(1L),
        "near_threshold_count"
      )))
    )
  )
}


.story_local_wcb_power_diagnostics <- function(
    threshold = 0.70,
    true_tau = 2.0,
    n_simulations = 50L,
    obs_per_cluster = 20L,
    seed = 42L,
    impose_null = TRUE,
    replay_seeds = NULL,
    weight_type = "rademacher"
) {
  .validate_story_local_wcb_power_diagnostics_args(
    threshold = threshold,
    true_tau = true_tau,
    n_simulations = n_simulations,
    obs_per_cluster = obs_per_cluster,
    seed = seed,
    impose_null = impose_null,
    replay_seeds = replay_seeds,
    weight_type = weight_type
  )

  threshold <- as.numeric(threshold)
  true_tau <- as.numeric(true_tau)
  n_simulations <- as.integer(n_simulations)
  obs_per_cluster <- as.integer(obs_per_cluster)
  seed <- as.integer(seed)
  weight_type <- as.character(weight_type)

  scenario <- list(
    id = "story_local_power_under_alternative",
    seed = seed,
    n_simulations = n_simulations,
    true_tau = true_tau,
    G = 10L,
    obs_per_cluster = obs_per_cluster,
    estimator = "wild_cluster_bootstrap",
    n_bootstrap = 199L,
    requested_n_bootstrap = 199L,
    weight_type = weight_type,
    alpha = 0.05,
    impose_null = isTRUE(impose_null),
    threshold = threshold,
    metric_name = "power_under_alternative",
    decision_direction = "lt",
    pvalue_threshold = 0.05
  )
  shared_dgp <- list(
    intercept = 10.0,
    cluster_effect_distribution = "normal(0, 2)",
    idiosyncratic_error_distribution = "normal(0, 1)",
    treatment_assignment = list(
      level = "cluster",
      default_probability = 0.5
    )
  )
  result <- .run_clustering_monte_carlo_wild_bootstrap_decision_scenario(
    scenario = scenario,
    shared_dgp = shared_dgp
  )
  result$rejection_count <- result$metric_hit_count
  result$power_confidence_interval <- result$metric_confidence_interval
  result$targeted_replays <- .build_story_local_wcb_targeted_replays(
    scenario = scenario,
    shared_dgp = shared_dgp,
    base_seed = seed,
    replay_seeds = replay_seeds,
    pvalue_threshold = result$pvalue_threshold
  )
  result$targeted_replay_summary <- .summarize_story_local_wcb_targeted_replays(
    result$targeted_replays
  )
  required_rejection_count <- as.integer(floor(threshold * n_simulations) + 1L)
  required_power_to_clear <- required_rejection_count / n_simulations
  rejection_count_shortfall <- max(
    required_rejection_count - result$rejection_count,
    0L
  )
  threshold_gap <- round(result$metric_value - threshold, digits = 12L)
  threshold_label <- formatC(threshold, format = "f", digits = 2)
  threshold_passed <- isTRUE(result$metric_value > threshold)
  blocker_boundary <- "story-local-hardening-power-gap-tc-9-4-19"
  runtime_blocker_boundary <- "d13-skip-backed-optional-backend-watchpoint"
  story_live_blockers <- if (threshold_passed) {
    runtime_blocker_boundary
  } else {
    c(blocker_boundary, runtime_blocker_boundary)
  }
  remaining_gap <- if (threshold_passed) {
    paste(
      sprintf(
        "TC-9.4.19 clears the %s story-local threshold;",
        threshold_label
      ),
      "dedicated-suite coverage is already landed, so the next step is to",
      "refresh the measured-blocker/state wording."
    )
  } else {
    paste(
      sprintf(
        "TC-9.4.19 currently measures power below the %s story-local threshold;",
        threshold_label
      ),
      "keep the dedicated-suite coverage landed and treat the case as an",
      "implementation hardening blocker."
    )
  }

  list(
    case_id = "TC-9.4.19",
    exact_status = if (threshold_passed) {
      "story-local-power-threshold-cleared"
    } else {
      "story-local-power-gap-measured"
    },
    blocker_boundary = blocker_boundary,
    story_live_blockers = story_live_blockers,
    d13_runtime_classification = runtime_blocker_boundary,
    global_blockers = if (threshold_passed) character(0) else blocker_boundary,
    dedicated_suite_ready = threshold_passed,
    remaining_gap = remaining_gap,
    scenario = scenario,
    result = result,
    threshold = threshold,
    required_rejection_count = required_rejection_count,
    required_power_to_clear = required_power_to_clear,
    rejection_count_shortfall = as.integer(rejection_count_shortfall),
    threshold_gap = threshold_gap,
    threshold_passed = threshold_passed
  )
}


.make_story_local_wcb_weight_family_result <- function(diagnostics) {
  list(
    weight_type = diagnostics$result$weight_type,
    rejection_count = as.integer(diagnostics$result$rejection_count),
    metric_value = as.numeric(diagnostics$result$metric_value),
    threshold_gap = as.numeric(diagnostics$threshold_gap),
    requested_n_bootstrap = as.integer(diagnostics$result$requested_n_bootstrap),
    actual_n_bootstrap = as.integer(diagnostics$result$actual_n_bootstrap),
    full_enumeration = isTRUE(diagnostics$result$full_enumeration),
    ci_method = diagnostics$result$ci_method,
    impose_null = isTRUE(diagnostics$result$impose_null),
    power_confidence_interval = diagnostics$result$power_confidence_interval,
    threshold_passed = isTRUE(diagnostics$threshold_passed)
  )
}


.story_local_wcb_weight_family_screening <- function() {
  weight_types <- c("rademacher", "mammen", "webb")
  diagnostics_by_weight <- lapply(weight_types, function(weight_type) {
    .story_local_wcb_power_diagnostics(weight_type = weight_type)
  })
  weight_results <- lapply(
    diagnostics_by_weight,
    .make_story_local_wcb_weight_family_result
  )
  metric_values <- vapply(weight_results, `[[`, numeric(1L), "metric_value")
  ranking <- vapply(weight_results, `[[`, character(1L), "weight_type")[
    order(-metric_values)
  ]
  best_non_rademacher <- weight_results[vapply(
    weight_results,
    function(entry) entry$weight_type != "rademacher",
    logical(1L)
  )]
  best_non_rademacher_index <- which.max(vapply(
    best_non_rademacher,
    `[[`,
    numeric(1L),
    "metric_value"
  ))
  best_non_rademacher_result <- best_non_rademacher[[best_non_rademacher_index]]
  rademacher_result <- weight_results[[1L]]
  rademacher_diagnostics <- diagnostics_by_weight[[1L]]

  list(
    case_id = rademacher_diagnostics$case_id,
    exact_status = "story-local-weight-family-screening-frozen",
    numeric_status = "simple-weight-switch-shortcuts-remain-below-threshold",
    blocker_boundary = rademacher_diagnostics$blocker_boundary,
    story_live_blockers = rademacher_diagnostics$story_live_blockers,
    consumer_summary_source = "story_local_wcb_power_diagnostics(weight_type = ...)",
    screened_weight_types = weight_types,
    weight_family_results = weight_results,
    comparison = list(
      ranking_best_to_worst = ranking,
      best_weight_type = ranking[[1L]],
      best_non_rademacher_weight_type = best_non_rademacher_result$weight_type,
      best_non_rademacher_minus_rademacher =
        best_non_rademacher_result$metric_value - rademacher_result$metric_value,
      webb_minus_mammen =
        weight_results[[3L]]$metric_value - weight_results[[2L]]$metric_value
    ),
    summary = list(
      canonical_rademacher_remains_best = identical(ranking[[1L]], "rademacher"),
      all_weight_families_fail_threshold = all(!vapply(
        weight_results,
        `[[`,
        logical(1L),
        "threshold_passed"
      )),
      non_rademacher_shortcut_eliminated =
        best_non_rademacher_result$metric_value < rademacher_result$metric_value,
      full_enumeration_only_for_rademacher =
        isTRUE(rademacher_result$full_enumeration) &&
        all(!vapply(weight_results[-1L], `[[`, logical(1L), "full_enumeration"))
    )
  )
}


.make_story_local_wcb_targeted_replay_case <- function(
    replication_row, base_seed, pvalue_threshold
) {
  list(
    attempt_id = as.integer(replication_row$attempt_id),
    replication_id = as.integer(replication_row$replication_id),
    seed = as.integer(base_seed + replication_row$attempt_id - 1L),
    replay_seed_offset = as.integer(replication_row$attempt_id - 1L),
    treated_cluster_count = as.integer(replication_row$treated_cluster_count),
    reject = isTRUE(replication_row$reject),
    pvalue = as.numeric(replication_row$pvalue),
    abs_gap_to_threshold =
      abs(as.numeric(replication_row$pvalue) - pvalue_threshold),
    att = as.numeric(replication_row$att),
    original_se = as.numeric(replication_row$original_se),
    t_stat_original = as.numeric(replication_row$t_stat_original)
  )
}


.story_local_wcb_targeted_replay_selector <- function() {
  diagnostics <- .story_local_wcb_power_diagnostics()
  replications <- diagnostics$result$replication_trace$replications
  selector_limit <- min(length(replications), 5L)

  replay_order <- order(
    vapply(
      replications,
      function(row) abs(as.numeric(row$pvalue) - diagnostics$result$pvalue_threshold),
      numeric(1L)
    ),
    vapply(replications, `[[`, integer(1L), "attempt_id")
  )
  targeted_replays <- lapply(
    replications[replay_order[seq_len(selector_limit)]],
    .make_story_local_wcb_targeted_replay_case,
    base_seed = as.integer(diagnostics$scenario$seed),
    pvalue_threshold = as.numeric(diagnostics$result$pvalue_threshold)
  )

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-targeted-replay-selector-contract-frozen",
    numeric_status = "targeted-replay-attempt-mapping-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "story_local_wcb_power_diagnostics()",
    base_seed = as.integer(diagnostics$scenario$seed),
    pvalue_threshold = as.numeric(diagnostics$result$pvalue_threshold),
    requested_n_bootstrap = as.integer(diagnostics$scenario$requested_n_bootstrap),
    actual_n_bootstrap = as.integer(diagnostics$result$actual_n_bootstrap),
    replay_seed_labels = vapply(
      targeted_replays,
      `[[`,
      integer(1L),
      "seed",
      USE.NAMES = FALSE
    ),
    targeted_replays = targeted_replays,
    summary = .summarize_story_local_wcb_targeted_replays(targeted_replays)
  )
}


.story_local_wcb_standard_t_counterfactual <- function() {
  build_counterfactual_cases <- function(rows, base_seed) {
    if (length(rows) == 0L) {
      return(list())
    }

    lapply(rows, function(row) {
      list(
        attempt_id = as.integer(row$attempt_id),
        replication_id = as.integer(row$replication_id),
        replay_seed = as.integer(base_seed + row$attempt_id - 1L),
        replay_seed_offset = as.integer(row$attempt_id - 1L),
        treated_cluster_count = as.integer(row$treated_cluster_count),
        wcb_pvalue = as.numeric(row$pvalue),
        standard_t_pvalue = as.numeric(row$standard_t_pvalue)
      )
    })
  }

  summarize_buckets <- function(rows) {
    if (length(rows) == 0L) {
      return(list())
    }

    bucket_counts <- table(vapply(
      rows,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ))
    ordered_names <- names(bucket_counts)[order(
      -as.integer(bucket_counts),
      as.integer(names(bucket_counts))
    )]

    lapply(ordered_names, function(bucket_name) {
      list(
        treated_cluster_count = as.integer(bucket_name),
        count = as.integer(bucket_counts[[bucket_name]])
      )
    })
  }

  diagnostics <- .story_local_wcb_power_diagnostics()
  replications <- diagnostics$result$replication_trace$replications
  counterfactual_summary <-
    diagnostics$result$replication_trace$standard_t_counterfactual_summary
  base_seed <- as.integer(diagnostics$scenario$seed)

  paired_rows <- replications[vapply(
    replications,
    function(row) {
      is.logical(row$standard_t_reject) &&
        length(row$standard_t_reject) == 1L &&
        !is.na(row$standard_t_reject) &&
        is.finite(row$standard_t_pvalue)
    },
    logical(1L)
  )]
  standard_t_reject <- vapply(
    paired_rows,
    `[[`,
    logical(1L),
    "standard_t_reject"
  )
  wild_bootstrap_reject <- vapply(
    paired_rows,
    `[[`,
    logical(1L),
    "reject"
  )
  standard_only_rows <- paired_rows[standard_t_reject & !wild_bootstrap_reject]
  wcb_only_rows <- paired_rows[!standard_t_reject & wild_bootstrap_reject]
  expected_standard_only_cases <- build_counterfactual_cases(
    standard_only_rows,
    base_seed = base_seed
  )
  expected_wcb_only_cases <- build_counterfactual_cases(
    wcb_only_rows,
    base_seed = base_seed
  )
  expected_standard_only_bucket_ranking <- summarize_buckets(standard_only_rows)
  expected_standard_only_seed_vector <- if (length(expected_standard_only_cases) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      expected_standard_only_cases,
      `[[`,
      integer(1L),
      "replay_seed"
    ))
  }
  summary_seed_vector <- as.integer(
    counterfactual_summary$standard_only_replay_seeds
  )
  standard_only_case_wcb_pvalues <- if (length(expected_standard_only_cases) == 0L) {
    numeric(0L)
  } else {
    vapply(
      expected_standard_only_cases,
      `[[`,
      numeric(1L),
      "wcb_pvalue"
    )
  }
  standard_only_case_standard_t_pvalues <- if (
    length(expected_standard_only_cases) == 0L
  ) {
    numeric(0L)
  } else {
    vapply(
      expected_standard_only_cases,
      `[[`,
      numeric(1L),
      "standard_t_pvalue"
    )
  }
  bucket_counts <- if (length(expected_standard_only_bucket_ranking) == 0L) {
    integer(0L)
  } else {
    vapply(
      expected_standard_only_bucket_ranking,
      `[[`,
      integer(1L),
      "count"
    )
  }
  bucket_labels <- if (length(expected_standard_only_bucket_ranking) == 0L) {
    integer(0L)
  } else {
    vapply(
      expected_standard_only_bucket_ranking,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    )
  }
  modal_standard_only_bucket <- if (length(bucket_counts) == 0L) {
    NA_integer_
  } else {
    as.integer(bucket_labels[[which.max(bucket_counts)]])
  }

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-standard-t-counterfactual-frozen",
    numeric_status = "paired-standard-t-vs-wcb-rescue-set-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "standard_t_counterfactual_summary",
    package_summary = counterfactual_summary,
    summary_consistency_checks = list(
      paired_count_matches_replications = identical(
        as.integer(counterfactual_summary$paired_replication_count),
        as.integer(length(paired_rows))
      ),
      standard_only_seed_vector_matches_cases = identical(
        summary_seed_vector,
        expected_standard_only_seed_vector
      ),
      standard_only_bucket_ranking_matches_cases = identical(
        counterfactual_summary$standard_only_bucket_ranking,
        expected_standard_only_bucket_ranking
      ),
      standard_only_cases_match_replications = identical(
        counterfactual_summary$standard_only_cases,
        expected_standard_only_cases
      ),
      wcb_only_cases_remain_empty = identical(
        counterfactual_summary$wcb_only_cases,
        expected_wcb_only_cases
      )
    ),
    paired_replication_count = as.integer(counterfactual_summary$paired_replication_count),
    standard_t_rejection_count =
      as.integer(counterfactual_summary$standard_t_rejection_count),
    wild_bootstrap_rejection_count =
      as.integer(counterfactual_summary$wild_bootstrap_rejection_count),
    shared_rejection_count = as.integer(counterfactual_summary$shared_rejection_count),
    standard_only_rejection_count =
      as.integer(counterfactual_summary$standard_only_rejection_count),
    wcb_only_rejection_count = as.integer(counterfactual_summary$wcb_only_rejection_count),
    paired_rejection_gap = as.integer(counterfactual_summary$paired_rejection_gap),
    standard_only_replay_seeds = as.integer(
      counterfactual_summary$standard_only_replay_seeds
    ),
    standard_only_near_threshold_replay_seeds = as.integer(
      counterfactual_summary$standard_only_near_threshold_replay_seeds
    ),
    standard_only_non_near_threshold_replay_seeds = as.integer(
      counterfactual_summary$standard_only_non_near_threshold_replay_seeds
    ),
    standard_only_bucket_ranking = counterfactual_summary$standard_only_bucket_ranking,
    standard_only_near_threshold_bucket_ranking =
      counterfactual_summary$standard_only_near_threshold_bucket_ranking,
    standard_only_non_near_threshold_bucket_ranking =
      counterfactual_summary$standard_only_non_near_threshold_bucket_ranking,
    standard_only_cases = counterfactual_summary$standard_only_cases,
    wcb_only_cases = counterfactual_summary$wcb_only_cases,
    standard_only_repair_frontier =
      counterfactual_summary$standard_only_repair_frontier,
    minimum_case_rank_to_clear_threshold =
      counterfactual_summary$minimum_case_rank_to_clear_threshold,
    summary = list(
      standard_t_power = as.numeric(
        counterfactual_summary$standard_t_rejection_count /
          counterfactual_summary$paired_replication_count
      ),
      wild_bootstrap_power = as.numeric(
        counterfactual_summary$wild_bootstrap_rejection_count /
          counterfactual_summary$paired_replication_count
      ),
      rejection_gap_share_of_replications = as.numeric(
        counterfactual_summary$paired_rejection_gap /
          counterfactual_summary$paired_replication_count
      ),
      standard_only_share_of_standard_rejections = as.numeric(
        counterfactual_summary$standard_only_rejection_count /
          counterfactual_summary$standard_t_rejection_count
      ),
      standard_only_share_of_wcb_non_rejections = as.numeric(
        counterfactual_summary$standard_only_rejection_count /
          (
            counterfactual_summary$paired_replication_count -
              counterfactual_summary$wild_bootstrap_rejection_count
          )
      ),
      standard_only_bucket_mode = modal_standard_only_bucket,
      standard_only_bucket_mode_count = if (length(bucket_counts) == 0L) {
        0L
      } else {
        as.integer(max(bucket_counts))
      },
      standard_only_near_threshold_share = as.numeric(
        counterfactual_summary$standard_only_near_threshold_share
      ),
      standard_only_non_near_threshold_share = as.numeric(
        counterfactual_summary$standard_only_non_near_threshold_share
      ),
      closest_standard_only_seed_to_threshold = as.integer(
        counterfactual_summary$closest_standard_only_seed_to_threshold
      ),
      closest_non_near_standard_only_seed = as.integer(
        counterfactual_summary$closest_non_near_standard_only_seed
      ),
      all_standard_only_cases_are_wcb_misses = all(
        standard_only_case_wcb_pvalues >= diagnostics$result$pvalue_threshold
      ),
      all_standard_only_cases_have_finite_standard_t_pvalues = all(
        is.finite(standard_only_case_standard_t_pvalues)
      )
    )
  )
}


.story_local_wcb_standard_t_rescue_threshold_split <- function() {
  annotate_cases <- function(cases, pvalue_threshold, near_threshold_upper) {
    if (length(cases) == 0L) {
      return(list())
    }

    lapply(cases, function(case) {
      c(
        case,
        list(
          wcb_gap_above_threshold = as.numeric(case$wcb_pvalue - pvalue_threshold),
          wcb_gap_above_window_upper = as.numeric(
            case$wcb_pvalue - near_threshold_upper
          ),
          standard_t_margin_below_threshold = as.numeric(
            pvalue_threshold - case$standard_t_pvalue
          )
        )
      )
    })
  }

  seed_vector_matches_cases <- function(cases, expected_seed_vector) {
    case_seed_vector <- if (length(cases) == 0L) {
      integer(0L)
    } else {
      as.integer(vapply(cases, `[[`, integer(1L), "replay_seed"))
    }

    identical(case_seed_vector, as.integer(expected_seed_vector))
  }

  diagnostics <- .story_local_wcb_power_diagnostics()
  counterfactual <- .story_local_wcb_standard_t_counterfactual()
  counterfactual_summary <- counterfactual$package_summary
  near_threshold_window <- counterfactual_summary$near_threshold_window
  pvalue_threshold <- as.numeric(diagnostics$result$pvalue_threshold)
  near_threshold_upper <- as.numeric(near_threshold_window$upper)
  standard_only_cases <- counterfactual_summary$standard_only_cases
  near_threshold_seed_vector <- as.integer(
    counterfactual_summary$standard_only_near_threshold_replay_seeds
  )
  non_near_threshold_seed_vector <- as.integer(
    counterfactual_summary$standard_only_non_near_threshold_replay_seeds
  )
  annotated_standard_only_cases <- annotate_cases(
    standard_only_cases,
    pvalue_threshold = pvalue_threshold,
    near_threshold_upper = near_threshold_upper
  )
  near_threshold_cases <- Filter(function(case) {
    case$replay_seed %in% near_threshold_seed_vector
  }, annotated_standard_only_cases)
  non_near_threshold_cases <- Filter(function(case) {
    case$replay_seed %in% non_near_threshold_seed_vector
  }, annotated_standard_only_cases)
  closest_standard_only_case <- if (length(annotated_standard_only_cases) == 0L) {
    NULL
  } else {
    annotated_standard_only_cases[[which.min(vapply(
      annotated_standard_only_cases,
      `[[`,
      numeric(1L),
      "wcb_gap_above_threshold"
    ))]]
  }
  closest_non_near_case <- if (length(non_near_threshold_cases) == 0L) {
    NULL
  } else {
    non_near_threshold_cases[[which.min(vapply(
      non_near_threshold_cases,
      `[[`,
      numeric(1L),
      "wcb_gap_above_window_upper"
    ))]]
  }
  paired_replication_count <- as.integer(counterfactual_summary$paired_replication_count)
  standard_only_rejection_count <- as.integer(
    counterfactual_summary$standard_only_rejection_count
  )
  standard_only_near_threshold_count <- as.integer(
    counterfactual_summary$standard_only_near_threshold_count
  )
  standard_only_non_near_threshold_count <- as.integer(
    counterfactual_summary$standard_only_non_near_threshold_count
  )

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-standard-t-rescue-threshold-split-frozen",
    numeric_status = "standard-only-rescues-threshold-split-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "standard_t_counterfactual_summary",
    package_summary = counterfactual_summary,
    near_threshold_window = near_threshold_window,
    standard_only_rejection_count = standard_only_rejection_count,
    standard_only_near_threshold_count = standard_only_near_threshold_count,
    standard_only_non_near_threshold_count =
      standard_only_non_near_threshold_count,
    standard_only_near_threshold_replay_seeds = near_threshold_seed_vector,
    standard_only_non_near_threshold_replay_seeds =
      non_near_threshold_seed_vector,
    standard_only_near_threshold_cases = near_threshold_cases,
    standard_only_non_near_threshold_cases = non_near_threshold_cases,
    standard_only_near_threshold_bucket_ranking =
      counterfactual_summary$standard_only_near_threshold_bucket_ranking,
    standard_only_non_near_threshold_bucket_ranking =
      counterfactual_summary$standard_only_non_near_threshold_bucket_ranking,
    standard_only_repair_frontier =
      counterfactual_summary$standard_only_repair_frontier,
    minimum_case_rank_to_clear_threshold =
      counterfactual_summary$minimum_case_rank_to_clear_threshold,
    summary_consistency_checks = list(
      standard_only_partition_is_exhaustive = identical(
        as.integer(length(near_threshold_cases) + length(non_near_threshold_cases)),
        as.integer(length(standard_only_cases))
      ),
      near_threshold_seed_vector_matches_cases = seed_vector_matches_cases(
        near_threshold_cases,
        near_threshold_seed_vector
      ),
      non_near_threshold_seed_vector_matches_cases = seed_vector_matches_cases(
        non_near_threshold_cases,
        non_near_threshold_seed_vector
      ),
      near_threshold_cases_match_window = all(vapply(
        near_threshold_cases,
        function(case) {
          case$wcb_pvalue >= as.numeric(near_threshold_window$lower) &&
            case$wcb_pvalue <= as.numeric(near_threshold_window$upper)
        },
        logical(1L)
      )),
      non_near_threshold_cases_exceed_window = all(vapply(
        non_near_threshold_cases,
        function(case) case$wcb_pvalue > as.numeric(near_threshold_window$upper),
        logical(1L)
      ))
    ),
    summary = list(
      near_threshold_share_of_standard_only = as.numeric(
        standard_only_near_threshold_count / standard_only_rejection_count
      ),
      non_near_threshold_share_of_standard_only = as.numeric(
        standard_only_non_near_threshold_count / standard_only_rejection_count
      ),
      near_threshold_share_of_replications = as.numeric(
        standard_only_near_threshold_count / paired_replication_count
      ),
      non_near_threshold_share_of_replications = as.numeric(
        standard_only_non_near_threshold_count / paired_replication_count
      ),
      closest_standard_only_seed_to_threshold = if (is.null(
        closest_standard_only_case
      )) {
        NULL
      } else {
        as.integer(closest_standard_only_case$replay_seed)
      },
      closest_standard_only_gap_to_threshold = if (is.null(
        closest_standard_only_case
      )) {
        NULL
      } else {
        as.numeric(closest_standard_only_case$wcb_gap_above_threshold)
      },
      closest_non_near_standard_only_seed = if (is.null(closest_non_near_case)) {
        NULL
      } else {
        as.integer(closest_non_near_case$replay_seed)
      },
      closest_non_near_gap_to_window_upper = if (is.null(closest_non_near_case)) {
        NULL
      } else {
        as.numeric(closest_non_near_case$wcb_gap_above_window_upper)
      }
    )
  )
}


.story_local_wcb_near_threshold_shortfall <- function() {
  diagnostics <- .story_local_wcb_power_diagnostics()
  replications <- diagnostics$result$replication_trace$replications
  near_threshold_window <- c(
    lower = diagnostics$result$pvalue_threshold - 0.01,
    upper = diagnostics$result$pvalue_threshold + 0.01
  )
  near_threshold_rows <- Filter(
    function(row) {
      row$pvalue >= near_threshold_window[["lower"]] &&
        row$pvalue <= near_threshold_window[["upper"]]
    },
    replications
  )
  near_threshold_total <- as.integer(length(near_threshold_rows))
  near_threshold_hits <- as.integer(sum(vapply(
    near_threshold_rows,
    function(row) isTRUE(row$reject),
    logical(1L)
  )))
  near_threshold_misses <- as.integer(near_threshold_total - near_threshold_hits)
  counterfactual_rejection_count <- as.integer(
    diagnostics$result$rejection_count + near_threshold_misses
  )
  counterfactual_power <- counterfactual_rejection_count / diagnostics$result$n_simulations
  residual_shortfall <- as.integer(max(
    diagnostics$required_rejection_count - counterfactual_rejection_count,
    0L
  ))

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-near-threshold-shortfall-bridge-frozen",
    numeric_status = "near-threshold-flip-capacity-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "story_local_wcb_power_diagnostics()",
    pvalue_threshold = as.numeric(diagnostics$result$pvalue_threshold),
    near_threshold_window = unname(as.numeric(near_threshold_window)),
    actual_rejection_count = as.integer(diagnostics$result$rejection_count),
    required_rejection_count = as.integer(diagnostics$required_rejection_count),
    rejection_count_shortfall = as.integer(diagnostics$rejection_count_shortfall),
    near_threshold_cases = list(
      total = near_threshold_total,
      hits = near_threshold_hits,
      misses = near_threshold_misses
    ),
    counterfactual_if_all_near_threshold_misses_flip = counterfactual_rejection_count,
    counterfactual_power_if_all_near_threshold_misses_flip =
      as.numeric(counterfactual_power),
    residual_shortfall_after_flipping_near_threshold_misses = residual_shortfall,
    can_clear_threshold_via_near_threshold_flips_alone =
      counterfactual_rejection_count >= diagnostics$required_rejection_count
  )
}


.story_local_wcb_standard_t_repair_frontier_replays <- function() {
  is_close_numeric <- function(lhs, rhs, tolerance = 1e-12) {
    if (length(lhs) != length(rhs)) {
      return(FALSE)
    }
    if (length(lhs) == 0L) {
      return(TRUE)
    }

    isTRUE(max(abs(lhs - rhs)) <= tolerance)
  }

  diagnostics <- .story_local_wcb_power_diagnostics()
  counterfactual_summary <-
    diagnostics$result$replication_trace$standard_t_counterfactual_summary
  frontier <- counterfactual_summary$standard_only_repair_frontier
  frontier_replay_seeds <- if (length(frontier) == 0L) {
    integer(0L)
  } else {
    vapply(frontier, `[[`, integer(1L), "replay_seed")
  }
  replay_diagnostics <- .story_local_wcb_power_diagnostics(
    replay_seeds = frontier_replay_seeds
  )
  targeted_replays <- replay_diagnostics$result$targeted_replays
  frontier_attempt_ids <- if (length(frontier) == 0L) {
    integer(0L)
  } else {
    vapply(frontier, `[[`, integer(1L), "attempt_id")
  }
  frontier_treated_cluster_pattern <- if (length(frontier) == 0L) {
    integer(0L)
  } else {
    vapply(frontier, `[[`, integer(1L), "treated_cluster_count")
  }
  frontier_wcb_pvalues <- if (length(frontier) == 0L) {
    numeric(0L)
  } else {
    vapply(frontier, `[[`, numeric(1L), "wcb_pvalue")
  }
  frontier_wcb_gaps <- if (length(frontier) == 0L) {
    numeric(0L)
  } else {
    vapply(frontier, `[[`, numeric(1L), "wcb_pvalue_gap")
  }
  targeted_replay_seeds <- if (length(targeted_replays) == 0L) {
    integer(0L)
  } else {
    vapply(targeted_replays, `[[`, integer(1L), "seed")
  }
  targeted_attempt_ids <- if (length(targeted_replays) == 0L) {
    integer(0L)
  } else {
    vapply(targeted_replays, `[[`, integer(1L), "attempt_id")
  }
  targeted_treated_cluster_pattern <- if (length(targeted_replays) == 0L) {
    integer(0L)
  } else {
    vapply(targeted_replays, `[[`, integer(1L), "treated_cluster_count")
  }
  targeted_replay_reject_pattern <- if (length(targeted_replays) == 0L) {
    logical(0L)
  } else {
    vapply(targeted_replays, `[[`, logical(1L), "reject")
  }
  targeted_replay_pvalues <- if (length(targeted_replays) == 0L) {
    numeric(0L)
  } else {
    vapply(targeted_replays, `[[`, numeric(1L), "pvalue")
  }
  targeted_replay_gaps <- if (length(targeted_replays) == 0L) {
    numeric(0L)
  } else {
    vapply(targeted_replays, `[[`, numeric(1L), "abs_gap_to_threshold")
  }

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-standard-t-repair-frontier-replays-frozen",
    numeric_status = "repair-frontier-targeted-replays-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "standard_t_counterfactual_summary",
    replay_selector_source =
      "story_local_wcb_power_diagnostics(replay_seeds = standard_only_repair_frontier)",
    package_summary = counterfactual_summary,
    standard_only_repair_frontier = frontier,
    targeted_replays = targeted_replays,
    summary_consistency_checks = list(
      targeted_replay_count_matches_frontier = identical(
        length(targeted_replays),
        length(frontier)
      ),
      replay_seed_order_matches_frontier = identical(
        targeted_replay_seeds,
        frontier_replay_seeds
      ),
      attempt_ids_match_frontier = identical(
        targeted_attempt_ids,
        frontier_attempt_ids
      ),
      treated_cluster_counts_match_frontier = identical(
        targeted_treated_cluster_pattern,
        frontier_treated_cluster_pattern
      ),
      all_targeted_replays_remain_non_rejections = all(
        !targeted_replay_reject_pattern
      ),
      wcb_pvalues_match_frontier = is_close_numeric(
        targeted_replay_pvalues,
        frontier_wcb_pvalues
      ),
      wcb_gap_matches_frontier = is_close_numeric(
        targeted_replay_gaps,
        frontier_wcb_gaps
      )
    ),
    summary = list(
      targeted_replay_count = as.integer(length(targeted_replays)),
      frontier_case_count = as.integer(length(frontier)),
      frontier_replay_seeds = frontier_replay_seeds,
      frontier_attempt_ids = frontier_attempt_ids,
      frontier_treated_cluster_pattern = frontier_treated_cluster_pattern,
      targeted_replay_reject_pattern = targeted_replay_reject_pattern,
      frontier_wcb_pvalues = frontier_wcb_pvalues,
      frontier_wcb_gaps = frontier_wcb_gaps,
      closest_frontier_seed = if (length(frontier_replay_seeds) == 0L) {
        NULL
      } else {
        as.integer(frontier_replay_seeds[[1L]])
      },
      furthest_frontier_seed = if (length(frontier_replay_seeds) == 0L) {
        NULL
      } else {
        as.integer(frontier_replay_seeds[[length(frontier_replay_seeds)]])
      }
    )
  )
}


.story_local_wcb_replication_mix_summary <- function() {
  diagnostics <- .story_local_wcb_power_diagnostics()
  profile <- diagnostics$result$replication_trace$treated_cluster_profile
  profile_summary <- diagnostics$result$replication_trace$treated_cluster_profile_summary
  modal_bucket <- profile_summary$modal_bucket
  contrast <- profile_summary$modal_bucket_signal_ratio_contrast
  anchor_bucket <- profile[[which(vapply(
    profile,
    function(bucket) {
      identical(
        bucket$treated_cluster_count,
        contrast$comparison_bucket_treated_cluster_count
      )
    },
    logical(1L)
  ))]]

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-power-gap-replication-mix-bridge-landed",
    numeric_status = "treated-cluster-mix-diagnostics-frozen",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "story_local_wcb_power_diagnostics()",
    threshold = as.numeric(diagnostics$threshold),
    replication_profile = profile,
    summary = list(
      total_replications = as.integer(diagnostics$result$n_simulations),
      total_rejections = as.integer(diagnostics$result$rejection_count),
      threshold_shortfall = as.integer(diagnostics$rejection_count_shortfall),
      modal_bucket = modal_bucket,
      signal_ratio_anchor_bucket = list(
        treated_cluster_count = as.integer(anchor_bucket$treated_cluster_count),
        avg_abs_att_to_se_ratio = as.numeric(anchor_bucket$avg_abs_att_to_se_ratio),
        median_abs_att_to_se_ratio = as.numeric(
          anchor_bucket$median_abs_att_to_se_ratio
        )
      ),
      signal_ratio_advantage_over_modal = list(
        avg_abs_att_to_se_ratio = as.numeric(
          anchor_bucket$avg_abs_att_to_se_ratio -
            modal_bucket$avg_abs_att_to_se_ratio
        ),
        median_abs_att_to_se_ratio = as.numeric(
          anchor_bucket$median_abs_att_to_se_ratio -
            modal_bucket$median_abs_att_to_se_ratio
        )
      ),
      high_rejection_buckets = profile_summary$high_rejection_buckets,
      zero_rejection_buckets = profile_summary$zero_rejection_buckets,
      near_threshold_buckets = profile_summary$near_threshold_buckets
    )
  )
}


.story_local_wcb_signal_ratio_bucket_profile <- function() {
  diagnostics <- .story_local_wcb_power_diagnostics()
  profile <- diagnostics$result$replication_trace$treated_cluster_profile
  profile_summary <- diagnostics$result$replication_trace$treated_cluster_profile_summary
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary
  contrast <- profile_summary$modal_bucket_signal_ratio_contrast
  profile_by_count <- stats::setNames(
    profile,
    vapply(
      profile,
      function(bucket) as.character(bucket$treated_cluster_count),
      character(1L)
    )
  )
  miss_bucket_ranking <- miss_summary$bucket_ranking
  miss_bucket_counts <- if (length(miss_bucket_ranking) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      miss_bucket_ranking,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ))
  }
  miss_counts_by_bucket <- stats::setNames(
    if (length(miss_bucket_ranking) == 0L) {
      integer(0L)
    } else {
      as.integer(vapply(
        miss_bucket_ranking,
        `[[`,
        integer(1L),
        "miss_count"
      ))
    },
    as.character(miss_bucket_counts)
  )
  miss_shares_by_bucket <- stats::setNames(
    if (length(miss_bucket_ranking) == 0L) {
      numeric(0L)
    } else {
      as.numeric(vapply(
        miss_bucket_ranking,
        `[[`,
        numeric(1L),
        "share_of_non_near_threshold_misses"
      ))
    },
    as.character(miss_bucket_counts)
  )
  miss_rank_by_bucket <- stats::setNames(
    as.integer(seq_along(miss_bucket_counts)),
    as.character(miss_bucket_counts)
  )
  signal_ratio_order <- order(
    -vapply(profile, `[[`, numeric(1L), "avg_abs_att_to_se_ratio"),
    vapply(profile, `[[`, integer(1L), "treated_cluster_count")
  )
  signal_ratio_ranking <- lapply(profile[signal_ratio_order], function(bucket) {
    bucket_key <- as.character(bucket$treated_cluster_count)

    list(
      treated_cluster_count = as.integer(bucket$treated_cluster_count),
      avg_abs_att_to_se_ratio = as.numeric(bucket$avg_abs_att_to_se_ratio),
      median_abs_att_to_se_ratio = as.numeric(bucket$median_abs_att_to_se_ratio),
      rejection_count = as.integer(bucket$rejection_count),
      rejection_rate = as.numeric(bucket$rejection_rate),
      avg_pvalue = as.numeric(bucket$avg_pvalue),
      near_threshold_count = as.integer(bucket$near_threshold_count),
      non_near_threshold_miss_count = if (bucket_key %in% names(miss_counts_by_bucket)) {
        as.integer(miss_counts_by_bucket[[bucket_key]])
      } else {
        0L
      },
      share_of_non_near_threshold_misses =
        if (bucket_key %in% names(miss_shares_by_bucket)) {
          as.numeric(miss_shares_by_bucket[[bucket_key]])
        } else {
          0
        },
      non_near_threshold_miss_rank = if (bucket_key %in% names(miss_rank_by_bucket)) {
        as.integer(miss_rank_by_bucket[[bucket_key]])
      } else {
        NA_integer_
      }
    )
  })
  signal_ratio_rank_counts <- as.integer(vapply(
    signal_ratio_ranking,
    `[[`,
    integer(1L),
    "treated_cluster_count"
  ))
  zero_rejection_bucket_counts <- as.integer(
    profile_summary$zero_rejection_buckets$treated_cluster_counts
  )
  high_rejection_bucket_counts <- as.integer(
    profile_summary$high_rejection_buckets$treated_cluster_counts
  )
  zero_rejection_buckets <- unname(
    profile_by_count[as.character(zero_rejection_bucket_counts)]
  )
  high_rejection_buckets <- unname(
    profile_by_count[as.character(high_rejection_bucket_counts)]
  )
  strongest_zero_rejection_bucket <- zero_rejection_buckets[[which.max(vapply(
    zero_rejection_buckets,
    `[[`,
    numeric(1L),
    "avg_abs_att_to_se_ratio"
  ))]]
  strongest_high_rejection_bucket <- high_rejection_buckets[[which.max(vapply(
    high_rejection_buckets,
    `[[`,
    numeric(1L),
    "avg_abs_att_to_se_ratio"
  ))]]
  strongest_bucket_overall_count <- signal_ratio_rank_counts[[1L]]
  strongest_bucket_overall_key <- as.character(strongest_bucket_overall_count)
  dominant_non_near_threshold_miss_bucket <- if (is.null(miss_summary$dominant_bucket)) {
    NULL
  } else {
    profile_by_count[[as.character(miss_summary$dominant_bucket$treated_cluster_count)]]
  }

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-signal-ratio-bucket-profile-frozen",
    numeric_status = "signal-ratio-bucket-profile-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "treated_cluster_profile_summary",
    profile_source = "treated_cluster_profile",
    signal_ratio_contrast = contrast,
    signal_ratio_ranking = signal_ratio_ranking,
    summary = list(
      modal_bucket = profile_summary$modal_bucket,
      comparison_bucket = profile_by_count[[
        as.character(contrast$comparison_bucket_treated_cluster_count)
      ]],
      strongest_bucket_overall = profile[[signal_ratio_order[[1L]]]],
      strongest_high_rejection_bucket = strongest_high_rejection_bucket,
      strongest_zero_rejection_bucket = strongest_zero_rejection_bucket,
      modal_bucket_rank_by_avg_ratio = as.integer(match(
        profile_summary$modal_bucket$treated_cluster_count,
        signal_ratio_rank_counts
      )),
      comparison_bucket_rank_by_avg_ratio = as.integer(match(
        contrast$comparison_bucket_treated_cluster_count,
        signal_ratio_rank_counts
      )),
      modal_vs_comparison_avg_ratio_gap = as.numeric(
        contrast$avg_ratio_gap_modal_minus_comparison
      ),
      modal_vs_comparison_median_ratio_gap = as.numeric(
        contrast$median_ratio_gap_modal_minus_comparison
      ),
      dominant_non_near_threshold_miss_bucket =
        dominant_non_near_threshold_miss_bucket,
      dominant_non_near_threshold_miss_bucket_rank_by_avg_ratio = if (
        !is.null(dominant_non_near_threshold_miss_bucket)
      ) {
        as.integer(match(
          dominant_non_near_threshold_miss_bucket$treated_cluster_count,
          signal_ratio_rank_counts
        ))
      } else {
        NA_integer_
      },
      dominant_non_near_threshold_miss_share = if (is.null(miss_summary$dominant_bucket)) {
        NaN
      } else {
        as.numeric(miss_summary$dominant_bucket$share_of_non_near_threshold_misses)
      },
      strongest_bucket_overall_non_near_threshold_miss_rank = if (
        strongest_bucket_overall_key %in% names(miss_rank_by_bucket)
      ) {
        as.integer(miss_rank_by_bucket[[strongest_bucket_overall_key]])
      } else {
        NA_integer_
      },
      strongest_bucket_overall_non_near_threshold_miss_share = if (
        strongest_bucket_overall_key %in% names(miss_shares_by_bucket)
      ) {
        as.numeric(miss_shares_by_bucket[[strongest_bucket_overall_key]])
      } else {
        0
      },
      high_rejection_bucket_counts = high_rejection_bucket_counts,
      zero_rejection_bucket_counts = zero_rejection_bucket_counts
    )
  )
}


.story_local_wcb_signal_ratio_repair_crosswalk <- function() {
  signal_ratio_profile <- .story_local_wcb_signal_ratio_bucket_profile()
  diagnostics <- .story_local_wcb_power_diagnostics()
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary
  repair_frontier <- miss_summary$bucket_repair_frontier

  signal_ratio_by_bucket <- stats::setNames(
    signal_ratio_profile$signal_ratio_ranking,
    as.character(vapply(
      signal_ratio_profile$signal_ratio_ranking,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ))
  )
  signal_ratio_rank_by_bucket <- stats::setNames(
    as.integer(seq_along(signal_ratio_profile$signal_ratio_ranking)),
    as.character(vapply(
      signal_ratio_profile$signal_ratio_ranking,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ))
  )
  repair_frontier_by_bucket <- stats::setNames(
    repair_frontier,
    as.character(vapply(
      repair_frontier,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ))
  )

  build_crosswalk_bucket <- function(treated_cluster_count) {
    bucket_key <- as.character(treated_cluster_count)
    signal_ratio_bucket <- signal_ratio_by_bucket[[bucket_key]]
    repair_bucket <- repair_frontier_by_bucket[[bucket_key]]

    c(
      signal_ratio_bucket,
      list(
        signal_ratio_rank = as.integer(signal_ratio_rank_by_bucket[[bucket_key]]),
        repair_frontier_rank = as.integer(repair_bucket$bucket_rank),
        cumulative_miss_count = as.integer(repair_bucket$cumulative_miss_count),
        cumulative_share_of_non_near_threshold_misses =
          as.numeric(repair_bucket$cumulative_share_of_non_near_threshold_misses),
        counterfactual_rejection_count =
          as.integer(repair_bucket$counterfactual_rejection_count),
        counterfactual_power = as.numeric(repair_bucket$counterfactual_power),
        residual_shortfall = as.integer(repair_bucket$residual_shortfall),
        clears_threshold = isTRUE(repair_bucket$clears_threshold),
        in_minimum_threshold_clearing_combo = treated_cluster_count %in%
          miss_summary$minimum_bucket_combo_to_clear_threshold
      )
    )
  }

  minimum_combo_counts <- as.integer(
    miss_summary$minimum_bucket_combo_to_clear_threshold
  )
  minimum_combo <- lapply(minimum_combo_counts, build_crosswalk_bucket)
  strongest_ratio_bucket_count <- as.integer(
    signal_ratio_profile$summary$strongest_bucket_overall$treated_cluster_count
  )
  strongest_ratio_bucket <- build_crosswalk_bucket(strongest_ratio_bucket_count)
  first_threshold_clearing_bucket <- minimum_combo[[length(minimum_combo)]]

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-signal-ratio-repair-crosswalk-frozen",
    numeric_status = "signal-ratio-repair-crosswalk-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "story_local_wcb_signal_ratio_bucket_profile()+bucket_repair_frontier",
    summary = list(
      strongest_ratio_bucket = strongest_ratio_bucket,
      minimum_threshold_clearing_combo = minimum_combo,
      first_threshold_clearing_bucket = first_threshold_clearing_bucket,
      strongest_ratio_bucket_excluded_from_repair_combo =
        !(strongest_ratio_bucket_count %in% minimum_combo_counts),
      threshold_clearing_combo_uses_mid_ratio_buckets = all(
        vapply(minimum_combo, `[[`, integer(1L), "signal_ratio_rank") >= 5L
      )
    )
  )
}


.story_local_wcb_non_near_threshold_miss_concentration <- function() {
  diagnostics <- .story_local_wcb_power_diagnostics()
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-non-near-threshold-miss-concentration-frozen",
    numeric_status = "non-near-threshold-miss-mass-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "story_local_wcb_power_diagnostics()",
    total_misses = as.integer(miss_summary$total_misses),
    near_threshold_misses = as.integer(miss_summary$near_threshold_misses),
    non_near_threshold_misses = as.integer(miss_summary$non_near_threshold_misses),
    non_near_threshold_share_of_misses = as.numeric(
      miss_summary$non_near_threshold_share_of_misses
    ),
    near_threshold_window = unname(as.numeric(miss_summary$near_threshold_window)),
    bucket_ranking = miss_summary$bucket_ranking,
    dominant_bucket = miss_summary$dominant_bucket,
    bucket_repair_frontier = miss_summary$bucket_repair_frontier,
    minimum_bucket_rank_to_clear_threshold =
      miss_summary$minimum_bucket_rank_to_clear_threshold,
    minimum_bucket_combo_to_clear_threshold =
      miss_summary$minimum_bucket_combo_to_clear_threshold,
    minimum_bucket_combo_share_of_non_near_threshold_misses =
      miss_summary$minimum_bucket_combo_share_of_non_near_threshold_misses,
    case_repair_frontier = miss_summary$case_repair_frontier,
    minimum_case_rank_to_clear_threshold =
      miss_summary$minimum_case_rank_to_clear_threshold,
    minimum_case_replay_seed_combo_to_clear_threshold =
      miss_summary$minimum_case_replay_seed_combo_to_clear_threshold,
    minimum_case_combo_share_of_non_near_threshold_misses =
      miss_summary$minimum_case_combo_share_of_non_near_threshold_misses
  )
}


.story_local_wcb_non_near_threshold_case_frontier <- function() {
  miss_concentration <- .story_local_wcb_non_near_threshold_miss_concentration()
  case_frontier <- miss_concentration$case_repair_frontier
  minimum_case_rank <- miss_concentration$minimum_case_rank_to_clear_threshold
  minimum_case_combo <- as.integer(
    miss_concentration$minimum_case_replay_seed_combo_to_clear_threshold
  )
  bucket_level_combo <- as.integer(
    miss_concentration$minimum_bucket_combo_to_clear_threshold
  )
  bucket_support <- if (is.null(minimum_case_rank) || minimum_case_rank == 0L) {
    integer(0L)
  } else {
    sort(unique(vapply(
      case_frontier[seq_len(minimum_case_rank)],
      `[[`,
      integer(1L),
      "treated_cluster_count"
    )))
  }
  bucket_frequency <- if (is.null(minimum_case_rank) || minimum_case_rank == 0L) {
    stats::setNames(integer(0L), character(0L))
  } else {
    combo_bucket_counts <- table(vapply(
      case_frontier[seq_len(minimum_case_rank)],
      function(row) as.character(row$treated_cluster_count),
      character(1L)
    ))
    ordered_bucket_names <- names(combo_bucket_counts)[order(
      -as.integer(combo_bucket_counts),
      as.integer(names(combo_bucket_counts))
    )]
    stats::setNames(
      as.integer(combo_bucket_counts[ordered_bucket_names]),
      ordered_bucket_names
    )
  }
  furthest_case <- if (is.null(minimum_case_rank) || minimum_case_rank == 0L) {
    NULL
  } else {
    case_frontier[[minimum_case_rank]]
  }

  list(
    case_id = miss_concentration$case_id,
    exact_status = "story-local-non-near-threshold-case-frontier-frozen",
    numeric_status = "case-level-threshold-clearing-frontier-machine-readable",
    blocker_boundary = miss_concentration$blocker_boundary,
    story_live_blockers = miss_concentration$story_live_blockers,
    consumer_summary_source = "non_near_threshold_miss_summary",
    total_misses = as.integer(miss_concentration$total_misses),
    non_near_threshold_misses = as.integer(
      miss_concentration$non_near_threshold_misses
    ),
    case_repair_frontier = case_frontier,
    minimum_case_rank_to_clear_threshold = minimum_case_rank,
    minimum_case_replay_seed_combo_to_clear_threshold = minimum_case_combo,
    minimum_case_combo_share_of_non_near_threshold_misses = as.numeric(
      miss_concentration$minimum_case_combo_share_of_non_near_threshold_misses
    ),
    summary = list(
      minimum_case_combo_bucket_support = as.integer(bucket_support),
      minimum_case_combo_bucket_frequency = bucket_frequency,
      minimum_case_combo_spans_more_buckets_than_bucket_combo =
        length(bucket_support) > length(bucket_level_combo),
      minimum_case_combo_uses_bucket_outside_bucket_level_combo = any(
        !bucket_support %in% bucket_level_combo
      ),
      furthest_case_seed_in_minimum_combo = if (is.null(furthest_case)) {
        NULL
      } else {
        as.integer(furthest_case$replay_seed)
      },
      furthest_case_gap_above_threshold_in_minimum_combo = if (
        is.null(furthest_case)
      ) {
        NULL
      } else {
        as.numeric(furthest_case$pvalue_gap_above_threshold)
      }
    )
  )
}


.story_local_wcb_non_near_threshold_case_frontier_tail_replays <- function() {
  case_frontier <- .story_local_wcb_non_near_threshold_case_frontier()
  minimum_case_rank <- case_frontier$minimum_case_rank_to_clear_threshold
  tail_size <- if (is.null(minimum_case_rank)) {
    0L
  } else {
    min(as.integer(minimum_case_rank), 5L)
  }
  tail_case_frontier <- if (tail_size == 0L) {
    list()
  } else {
    case_frontier$case_repair_frontier[
      seq.int(minimum_case_rank - tail_size + 1L, minimum_case_rank)
    ]
  }
  tail_case_replay_seeds <- if (tail_size == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      tail_case_frontier,
      `[[`,
      integer(1L),
      "replay_seed"
    ))
  }
  replay_diagnostics <- .story_local_wcb_power_diagnostics(
    replay_seeds = tail_case_replay_seeds
  )
  targeted_replays <- replay_diagnostics$result$targeted_replays
  tail_case_frontier <- lapply(seq_along(tail_case_frontier), function(idx) {
    c(
      tail_case_frontier[[idx]],
      list(targeted_replay = targeted_replays[[idx]])
    )
  })
  pre_threshold_case <- if (length(tail_case_frontier) >= 2L) {
    tail_case_frontier[[length(tail_case_frontier) - 1L]]
  } else {
    NULL
  }
  threshold_crossing_case <- if (length(tail_case_frontier) >= 1L) {
    tail_case_frontier[[length(tail_case_frontier)]]
  } else {
    NULL
  }
  tail_gaps <- if (length(tail_case_frontier) == 0L) {
    numeric(0L)
  } else {
    vapply(
      tail_case_frontier,
      `[[`,
      numeric(1L),
      "pvalue_gap_above_threshold"
    )
  }
  lowest_gap_tail_seed <- if (length(tail_case_frontier) == 0L) {
    NULL
  } else {
    tail_case_frontier[[which.min(tail_gaps)]]$replay_seed
  }

  list(
    case_id = case_frontier$case_id,
    exact_status = "story-local-non-near-threshold-case-frontier-tail-replays-frozen",
    numeric_status = "threshold-edge-case-tail-replays-machine-readable",
    blocker_boundary = case_frontier$blocker_boundary,
    story_live_blockers = case_frontier$story_live_blockers,
    consumer_summary_source = "story_local_wcb_non_near_threshold_case_frontier()",
    replay_selector_source =
      "story_local_wcb_power_diagnostics(replay_seeds = threshold_edge_case_tail)",
    tail_case_replay_seeds = tail_case_replay_seeds,
    tail_case_frontier = tail_case_frontier,
    summary = list(
      all_tail_targeted_replays_remain_misses = all(!vapply(
        lapply(tail_case_frontier, `[[`, "targeted_replay"),
        `[[`,
        logical(1L),
        "reject"
      )),
      pre_threshold_case_rank = if (is.null(pre_threshold_case)) {
        NULL
      } else {
        as.integer(pre_threshold_case$case_rank)
      },
      pre_threshold_replay_seed = if (is.null(pre_threshold_case)) {
        NULL
      } else {
        as.integer(pre_threshold_case$replay_seed)
      },
      pre_threshold_counterfactual_power = if (is.null(pre_threshold_case)) {
        NULL
      } else {
        as.numeric(pre_threshold_case$counterfactual_power)
      },
      threshold_crossing_case_rank = if (is.null(threshold_crossing_case)) {
        NULL
      } else {
        as.integer(threshold_crossing_case$case_rank)
      },
      threshold_crossing_replay_seed = if (is.null(threshold_crossing_case)) {
        NULL
      } else {
        as.integer(threshold_crossing_case$replay_seed)
      },
      threshold_crossing_seed_still_observed_miss = if (
        is.null(threshold_crossing_case)
      ) {
        FALSE
      } else {
        !isTRUE(threshold_crossing_case$targeted_replay$reject)
      },
      lowest_gap_tail_seed = if (is.null(lowest_gap_tail_seed)) {
        NULL
      } else {
        as.integer(lowest_gap_tail_seed)
      },
      threshold_crossing_seed_has_largest_gap_in_tail = if (
        is.null(threshold_crossing_case) || length(tail_gaps) == 0L
      ) {
        FALSE
      } else {
        identical(
          as.integer(threshold_crossing_case$replay_seed),
          as.integer(tail_case_frontier[[which.max(tail_gaps)]]$replay_seed)
        )
      }
    )
  )
}


.story_local_wcb_non_near_threshold_case_frontier_residual_replays <- function() {
  case_frontier <- .story_local_wcb_non_near_threshold_case_frontier()
  minimum_case_rank <- case_frontier$minimum_case_rank_to_clear_threshold
  full_case_frontier <- case_frontier$case_repair_frontier
  residual_case_frontier <- if (
    is.null(minimum_case_rank) ||
      minimum_case_rank >= length(full_case_frontier)
  ) {
    list()
  } else {
    full_case_frontier[
      seq.int(minimum_case_rank + 1L, length(full_case_frontier))
    ]
  }
  residual_case_replay_seeds <- if (length(residual_case_frontier) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      residual_case_frontier,
      `[[`,
      integer(1L),
      "replay_seed"
    ))
  }
  replay_diagnostics <- .story_local_wcb_power_diagnostics(
    replay_seeds = residual_case_replay_seeds
  )
  targeted_replays <- replay_diagnostics$result$targeted_replays
  residual_case_frontier <- lapply(seq_along(residual_case_frontier), function(idx) {
    c(
      residual_case_frontier[[idx]],
      list(targeted_replay = targeted_replays[[idx]])
    )
  })
  residual_bucket_support <- if (length(residual_case_frontier) == 0L) {
    integer(0L)
  } else {
    sort(unique(vapply(
      residual_case_frontier,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    )))
  }
  residual_bucket_frequency <- if (length(residual_case_frontier) == 0L) {
    stats::setNames(integer(0L), character(0L))
  } else {
    bucket_counts <- table(vapply(
      residual_case_frontier,
      function(row) as.character(row$treated_cluster_count),
      character(1L)
    ))
    ordered_bucket_names <- names(bucket_counts)[order(
      -as.integer(bucket_counts),
      as.integer(names(bucket_counts))
    )]
    stats::setNames(
      as.integer(bucket_counts[ordered_bucket_names]),
      ordered_bucket_names
    )
  }
  residual_gaps <- if (length(residual_case_frontier) == 0L) {
    numeric(0L)
  } else {
    vapply(
      residual_case_frontier,
      `[[`,
      numeric(1L),
      "pvalue_gap_above_threshold"
    )
  }
  easiest_residual_case <- if (length(residual_case_frontier) == 0L) {
    NULL
  } else {
    residual_case_frontier[[which.min(residual_gaps)]]
  }
  hardest_residual_case <- if (length(residual_case_frontier) == 0L) {
    NULL
  } else {
    residual_case_frontier[[which.max(residual_gaps)]]
  }
  minimum_case_bucket_support <- as.integer(
    case_frontier$summary$minimum_case_combo_bucket_support
  )

  list(
    case_id = case_frontier$case_id,
    exact_status = "story-local-non-near-threshold-case-frontier-residual-replays-frozen",
    numeric_status = "post-threshold-residual-replays-machine-readable",
    blocker_boundary = case_frontier$blocker_boundary,
    story_live_blockers = case_frontier$story_live_blockers,
    consumer_summary_source = "story_local_wcb_non_near_threshold_case_frontier()",
    replay_selector_source =
      "story_local_wcb_power_diagnostics(replay_seeds = residual_case_frontier)",
    residual_case_replay_seeds = residual_case_replay_seeds,
    residual_case_frontier = residual_case_frontier,
    summary = list(
      residual_case_count = as.integer(length(residual_case_frontier)),
      threshold_clearing_rank = if (is.null(minimum_case_rank)) {
        NULL
      } else {
        as.integer(minimum_case_rank)
      },
      first_residual_case_rank = if (length(residual_case_frontier) == 0L) {
        NULL
      } else {
        as.integer(residual_case_frontier[[1L]]$case_rank)
      },
      easiest_residual_replay_seed = if (is.null(easiest_residual_case)) {
        NULL
      } else {
        as.integer(easiest_residual_case$replay_seed)
      },
      easiest_residual_gap_above_threshold = if (is.null(easiest_residual_case)) {
        NULL
      } else {
        as.numeric(easiest_residual_case$pvalue_gap_above_threshold)
      },
      hardest_residual_replay_seed = if (is.null(hardest_residual_case)) {
        NULL
      } else {
        as.integer(hardest_residual_case$replay_seed)
      },
      hardest_residual_gap_above_threshold = if (is.null(hardest_residual_case)) {
        NULL
      } else {
        as.numeric(hardest_residual_case$pvalue_gap_above_threshold)
      },
      residual_bucket_support = as.integer(residual_bucket_support),
      residual_bucket_frequency = residual_bucket_frequency,
      all_residual_targeted_replays_remain_misses = all(!vapply(
        lapply(residual_case_frontier, `[[`, "targeted_replay"),
        `[[`,
        logical(1L),
        "reject"
      )),
      bucket7_only_enters_after_threshold = (
        any(residual_bucket_support == 7L) &&
          !any(minimum_case_bucket_support == 7L)
      ),
      residual_family_excludes_bucket8_and_9 = !any(
        residual_bucket_support %in% c(8L, 9L)
      )
    )
  )
}


.story_local_wcb_case_combo_residual_partition <- function() {
  case_frontier <- .story_local_wcb_non_near_threshold_case_frontier()
  residual_replays <- .story_local_wcb_non_near_threshold_case_frontier_residual_replays()

  minimum_case_rank <- if (is.null(case_frontier$minimum_case_rank_to_clear_threshold)) {
    0L
  } else {
    as.integer(case_frontier$minimum_case_rank_to_clear_threshold)
  }
  non_near_threshold_misses <- as.integer(case_frontier$non_near_threshold_misses)
  minimum_bucket_support <- as.integer(case_frontier$summary$minimum_case_combo_bucket_support)
  residual_bucket_support <- as.integer(residual_replays$summary$residual_bucket_support)
  overlap_bucket_support <- intersect(minimum_bucket_support, residual_bucket_support)
  minimum_only_bucket_support <- setdiff(minimum_bucket_support, residual_bucket_support)
  residual_only_bucket_support <- setdiff(residual_bucket_support, minimum_bucket_support)
  residual_case_count <- as.integer(residual_replays$summary$residual_case_count)
  first_residual_case_rank <- residual_replays$summary$first_residual_case_rank

  list(
    case_id = case_frontier$case_id,
    exact_status = "story-local-case-combo-residual-partition-frozen",
    numeric_status = "threshold-clearing-vs-residual-frontier-machine-readable",
    blocker_boundary = case_frontier$blocker_boundary,
    story_live_blockers = case_frontier$story_live_blockers,
    consumer_summary_source = paste0(
      "story_local_wcb_non_near_threshold_case_frontier()+",
      "story_local_wcb_non_near_threshold_case_frontier_residual_replays()"
    ),
    minimum_case_rank_to_clear_threshold = minimum_case_rank,
    minimum_case_combo_replay_seeds = as.integer(
      case_frontier$minimum_case_replay_seed_combo_to_clear_threshold
    ),
    residual_case_replay_seeds = as.integer(residual_replays$residual_case_replay_seeds),
    overlap_bucket_support = as.integer(overlap_bucket_support),
    minimum_only_bucket_support = as.integer(minimum_only_bucket_support),
    residual_only_bucket_support = as.integer(residual_only_bucket_support),
    summary = list(
      non_near_threshold_miss_count = non_near_threshold_misses,
      minimum_case_combo_count = minimum_case_rank,
      residual_case_count = residual_case_count,
      minimum_case_combo_share_of_non_near_threshold_misses = if (
        non_near_threshold_misses > 0L
      ) {
        minimum_case_rank / non_near_threshold_misses
      } else {
        NaN
      },
      residual_case_share_of_non_near_threshold_misses = if (
        non_near_threshold_misses > 0L
      ) {
        residual_case_count / non_near_threshold_misses
      } else {
        NaN
      },
      combined_partition_reconstructs_non_near_threshold_misses = (
        minimum_case_rank + residual_case_count
      ) == non_near_threshold_misses,
      residual_frontier_starts_after_threshold = !is.null(first_residual_case_rank) &&
        identical(as.integer(first_residual_case_rank), minimum_case_rank + 1L),
      bucket7_is_residual_only = length(residual_only_bucket_support) == 1L &&
        identical(as.integer(residual_only_bucket_support[[1L]]), 7L),
      buckets8_and_9_are_threshold_only = identical(
        as.integer(minimum_only_bucket_support),
        c(8L, 9L)
      )
    )
  )
}


.story_local_wcb_power_repair_ladder <- function() {
  diagnostics <- .story_local_wcb_power_diagnostics()
  near_threshold_shortfall <- .story_local_wcb_near_threshold_shortfall()
  standard_t_counterfactual <- .story_local_wcb_standard_t_counterfactual()
  miss_concentration <- .story_local_wcb_non_near_threshold_miss_concentration()
  minimum_case_combo_seeds <- as.integer(
    miss_concentration$minimum_case_replay_seed_combo_to_clear_threshold
  )
  minimum_bucket_combo_counts <- as.integer(
    miss_concentration$minimum_bucket_combo_to_clear_threshold
  )
  case_repair_frontier <- miss_concentration$case_repair_frontier
  first_case_clear_rung <- if (
    is.null(miss_concentration$minimum_case_rank_to_clear_threshold)
  ) {
    NULL
  } else {
    case_repair_frontier[[miss_concentration$minimum_case_rank_to_clear_threshold]]
  }
  repair_frontier <- miss_concentration$bucket_repair_frontier
  first_clear_idx <- which(vapply(
    repair_frontier,
    function(bucket) isTRUE(bucket$clears_threshold),
    logical(1L)
  ))
  first_clear_rung <- if (length(first_clear_idx) == 0L) {
    NULL
  } else {
    repair_frontier[[first_clear_idx[[1L]]]]
  }

  repair_rungs <- list(
    list(
      rung_id = "observed-wild-bootstrap",
      rejection_count = as.integer(diagnostics$result$rejection_count),
      power = as.numeric(diagnostics$result$metric_value),
      residual_shortfall = as.integer(diagnostics$rejection_count_shortfall),
      clears_threshold = isTRUE(diagnostics$threshold_passed)
    ),
    list(
      rung_id = "all-near-threshold-misses-flip",
      rejection_count = as.integer(
        near_threshold_shortfall$counterfactual_if_all_near_threshold_misses_flip
      ),
      power = as.numeric(
        near_threshold_shortfall$counterfactual_power_if_all_near_threshold_misses_flip
      ),
      residual_shortfall = as.integer(
        near_threshold_shortfall$residual_shortfall_after_flipping_near_threshold_misses
      ),
      clears_threshold = isTRUE(
        near_threshold_shortfall$can_clear_threshold_via_near_threshold_flips_alone
      )
    ),
    list(
      rung_id = "standard-t-counterfactual",
      rejection_count = as.integer(
        standard_t_counterfactual$standard_t_rejection_count
      ),
      power = as.numeric(standard_t_counterfactual$summary$standard_t_power),
      residual_shortfall = as.integer(max(
        diagnostics$required_rejection_count -
          standard_t_counterfactual$standard_t_rejection_count,
        0L
      )),
      clears_threshold = isTRUE(
        standard_t_counterfactual$standard_t_rejection_count >=
          diagnostics$required_rejection_count
      )
    ),
    list(
      rung_id = "minimum-non-near-threshold-case-combo",
      rejection_count = if (is.null(first_case_clear_rung)) {
        NULL
      } else {
        as.integer(first_case_clear_rung$counterfactual_rejection_count)
      },
      power = if (is.null(first_case_clear_rung)) {
        NULL
      } else {
        as.numeric(first_case_clear_rung$counterfactual_power)
      },
      residual_shortfall = if (is.null(first_case_clear_rung)) {
        NULL
      } else {
        as.integer(first_case_clear_rung$residual_shortfall)
      },
      clears_threshold = if (is.null(first_case_clear_rung)) {
        FALSE
      } else {
        isTRUE(first_case_clear_rung$clears_threshold)
      },
      replay_seed_combo = minimum_case_combo_seeds,
      share_of_non_near_threshold_misses = as.numeric(
        miss_concentration$minimum_case_combo_share_of_non_near_threshold_misses
      )
    ),
    list(
      rung_id = "minimum-non-near-threshold-bucket-combo",
      rejection_count = if (is.null(first_clear_rung)) {
        NULL
      } else {
        as.integer(first_clear_rung$counterfactual_rejection_count)
      },
      power = if (is.null(first_clear_rung)) {
        NULL
      } else {
        as.numeric(first_clear_rung$counterfactual_power)
      },
      residual_shortfall = if (is.null(first_clear_rung)) {
        NULL
      } else {
        as.integer(first_clear_rung$residual_shortfall)
      },
      clears_threshold = if (is.null(first_clear_rung)) {
        FALSE
      } else {
        isTRUE(first_clear_rung$clears_threshold)
      },
      treated_cluster_counts = minimum_bucket_combo_counts,
      share_of_non_near_threshold_misses = as.numeric(
        miss_concentration$minimum_bucket_combo_share_of_non_near_threshold_misses
      )
    )
  )
  first_threshold_clearing_rung_idx <- which(vapply(
    repair_rungs,
    function(rung) isTRUE(rung$clears_threshold),
    logical(1L)
  ))

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-power-repair-ladder-frozen",
    numeric_status = "power-repair-ladder-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = paste0(
      "story_local_wcb_near_threshold_shortfall()+",
      "story_local_wcb_standard_t_counterfactual()+",
      "story_local_wcb_non_near_threshold_miss_concentration()"
    ),
    threshold = as.numeric(diagnostics$threshold),
    required_rejection_count = as.integer(diagnostics$required_rejection_count),
    required_power_to_clear = as.numeric(diagnostics$required_power_to_clear),
    repair_rungs = repair_rungs,
    summary = list(
      first_threshold_clearing_rung = if (
        length(first_threshold_clearing_rung_idx) == 0L
      ) {
        NULL
      } else {
        repair_rungs[[first_threshold_clearing_rung_idx[[1L]]]]$rung_id
      },
      minimum_threshold_clearing_case_combo = minimum_case_combo_seeds,
      minimum_threshold_clearing_bucket_combo = minimum_bucket_combo_counts,
      near_threshold_alone_still_below_standard_t_counterfactual = (
        near_threshold_shortfall$counterfactual_power_if_all_near_threshold_misses_flip <
          standard_t_counterfactual$summary$standard_t_power
      ),
      standard_t_counterfactual_still_below_threshold = (
        standard_t_counterfactual$summary$standard_t_power <
          diagnostics$required_power_to_clear
      ),
      minimum_case_combo_clears_with_fewer_repairs_than_bucket_combo = (
        !is.null(first_case_clear_rung) &&
          !is.null(first_clear_rung) &&
          length(minimum_case_combo_seeds) < first_clear_rung$cumulative_miss_count
      ),
      threshold_crossing_requires_non_near_threshold_combo = (
        !isTRUE(
          near_threshold_shortfall$can_clear_threshold_via_near_threshold_flips_alone
        ) &&
          isTRUE(repair_rungs[[4L]]$clears_threshold)
      )
    )
  )
}


.story_local_wcb_top_bucket_representative_replays <- function() {
  diagnostics <- .story_local_wcb_power_diagnostics()
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary
  top_bucket_limit <- min(length(miss_summary$bucket_ranking), 3L)
  top_buckets <- if (top_bucket_limit > 0L) {
    miss_summary$bucket_ranking[seq_len(top_bucket_limit)]
  } else {
    list()
  }
  representative_replay_seeds <- if (length(top_buckets) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      top_buckets,
      function(bucket) bucket$representative_large_miss$replay_seed,
      integer(1L)
    ))
  }
  replay_diagnostics <- .story_local_wcb_power_diagnostics(
    replay_seeds = representative_replay_seeds
  )
  targeted_replays <- replay_diagnostics$result$targeted_replays

  top_bucket_replays <- lapply(seq_along(top_buckets), function(idx) {
    bucket <- top_buckets[[idx]]

    list(
      treated_cluster_count = as.integer(bucket$treated_cluster_count),
      miss_count = as.integer(bucket$miss_count),
      share_of_non_near_threshold_misses = as.numeric(
        bucket$share_of_non_near_threshold_misses
      ),
      representative_large_miss = bucket$representative_large_miss,
      targeted_replay = targeted_replays[[idx]]
    )
  })
  top_bucket_miss_counts <- if (length(top_bucket_replays) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      top_bucket_replays,
      `[[`,
      integer(1L),
      "miss_count"
    ))
  }
  top_bucket_order <- if (length(top_bucket_replays) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      top_bucket_replays,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ))
  }

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-top-bucket-representative-replays-frozen",
    numeric_status = "top-bucket-far-miss-replays-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "non_near_threshold_miss_summary",
    replay_selector_source =
      "story_local_wcb_power_diagnostics(replay_seeds = representative_large_miss$replay_seed)",
    total_misses = as.integer(miss_summary$total_misses),
    non_near_threshold_misses = as.integer(
      miss_summary$non_near_threshold_misses
    ),
    top_bucket_order = top_bucket_order,
    representative_replay_seeds = representative_replay_seeds,
    top_bucket_replays = top_bucket_replays,
    summary = list(
      top_three_share_of_non_near_threshold_misses = if (
        miss_summary$non_near_threshold_misses > 0L
      ) {
        sum(top_bucket_miss_counts) /
          as.numeric(miss_summary$non_near_threshold_misses)
      } else {
        NaN
      },
      bucket5_share_of_top_three_misses = if (sum(top_bucket_miss_counts) > 0L) {
        top_bucket_miss_counts[[1L]] / sum(top_bucket_miss_counts)
      } else {
        NaN
      },
      bucket5_excess_miss_count_vs_bucket6 = if (
        length(top_bucket_miss_counts) >= 2L
      ) {
        as.integer(top_bucket_miss_counts[[1L]] - top_bucket_miss_counts[[2L]])
      } else {
        NULL
      },
      bucket5_excess_miss_count_vs_bucket4 = if (
        length(top_bucket_miss_counts) >= 3L
      ) {
        as.integer(top_bucket_miss_counts[[1L]] - top_bucket_miss_counts[[3L]])
      } else {
        NULL
      },
      all_replayed_cases_remain_misses = all(!vapply(
        targeted_replays,
        `[[`,
        logical(1L),
        "reject"
      )),
      targeted_replay_summary = replay_diagnostics$result$targeted_replay_summary
    )
  )
}


.story_local_wcb_representative_far_miss_mechanism_summary <- function() {
  diagnostics <- .story_local_wcb_power_diagnostics()
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary
  bucket_ranking <- miss_summary$bucket_ranking
  top_bucket_limit <- min(length(bucket_ranking), 3L)
  top_buckets <- if (top_bucket_limit > 0L) {
    bucket_ranking[seq_len(top_bucket_limit)]
  } else {
    list()
  }
  mechanism_summary <- miss_summary$representative_far_miss_mechanism_summary
  bucket_reference_by_count <- stats::setNames(
    top_buckets,
    vapply(
      top_buckets,
      function(bucket) as.character(bucket$treated_cluster_count),
      character(1L)
    )
  )

  bucket_mechanisms <- lapply(
    mechanism_summary$bucket_mechanisms,
    function(bucket_summary) {
      bucket <- bucket_reference_by_count[[
        as.character(bucket_summary$treated_cluster_count)
      ]]
      representative_row <- bucket$representative_large_miss

      list(
        treated_cluster_count = as.integer(bucket_summary$treated_cluster_count),
        representative_replay_seed = as.integer(
          bucket_summary$representative_replay_seed
        ),
        mechanism_label = bucket_summary$mechanism_label,
        avg_abs_att = as.numeric(bucket$avg_abs_att),
        avg_original_se = as.numeric(bucket$avg_original_se),
        avg_abs_t_stat = as.numeric(bucket$avg_abs_t_stat),
        representative_abs_att = as.numeric(representative_row$abs_att),
        representative_original_se = as.numeric(representative_row$original_se),
        representative_abs_att_to_bucket_avg_ratio = as.numeric(
          bucket_summary$representative_abs_att_to_bucket_avg_ratio
        ),
        representative_original_se_to_bucket_avg_ratio = as.numeric(
          bucket_summary$representative_original_se_to_bucket_avg_ratio
        ),
        representative_abs_t_stat_to_bucket_avg_ratio = as.numeric(
          representative_row$abs_t_stat_to_bucket_avg_ratio
        )
      )
    }
  )

  bucket_order_from_bucket_ranking <- if (length(top_buckets) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      top_buckets,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ))
  }
  replay_seeds_from_bucket_ranking <- if (length(top_buckets) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      top_buckets,
      function(bucket) bucket$representative_large_miss$replay_seed,
      integer(1L)
    ))
  }
  mechanism_labels_from_bucket_ranking <- if (length(top_buckets) == 0L) {
    character(0L)
  } else {
    vapply(
      top_buckets,
      function(bucket) {
        .classify_story_local_wcb_representative_far_miss_mechanism(
          bucket$representative_large_miss
        )
      },
      character(1L)
    )
  }
  ratio_matches_bucket_ranking <- all(vapply(
    seq_along(bucket_mechanisms),
    function(idx) {
      bucket <- top_buckets[[idx]]
      bucket_mechanism <- bucket_mechanisms[[idx]]

      isTRUE(all.equal(
        bucket_mechanism$representative_abs_att_to_bucket_avg_ratio,
        as.numeric(bucket$representative_large_miss$abs_att_to_bucket_avg_ratio),
        tolerance = 1e-12
      )) && isTRUE(all.equal(
        bucket_mechanism$representative_original_se_to_bucket_avg_ratio,
        as.numeric(bucket$representative_large_miss$original_se_to_bucket_avg_ratio),
        tolerance = 1e-12
      ))
    },
    logical(1L)
  ))

  list(
    case_id = diagnostics$case_id,
    exact_status = "story-local-representative-far-miss-mechanism-frozen",
    numeric_status = "representative-far-miss-signal-vs-se-mechanism-machine-readable",
    blocker_boundary = diagnostics$blocker_boundary,
    story_live_blockers = diagnostics$story_live_blockers,
    consumer_summary_source = "representative_far_miss_mechanism_summary",
    package_summary = mechanism_summary,
    summary_consistency_checks = list(
      top_bucket_order_matches_bucket_ranking = identical(
        as.integer(mechanism_summary$top_bucket_order),
        bucket_order_from_bucket_ranking
      ),
      representative_replay_seeds_match_bucket_ranking = identical(
        as.integer(mechanism_summary$representative_replay_seeds),
        replay_seeds_from_bucket_ranking
      ),
      mechanism_labels_match_bucket_ranking = identical(
        vapply(bucket_mechanisms, `[[`, character(1L), "mechanism_label"),
        mechanism_labels_from_bucket_ranking
      ),
      ratios_match_bucket_ranking = isTRUE(ratio_matches_bucket_ranking)
    ),
    top_bucket_order = as.integer(mechanism_summary$top_bucket_order),
    representative_replay_seeds = as.integer(
      mechanism_summary$representative_replay_seeds
    ),
    low_signal_plus_elevated_se_buckets = as.integer(
      mechanism_summary$low_signal_plus_elevated_se_buckets
    ),
    low_signal_dominant_buckets = as.integer(
      mechanism_summary$low_signal_dominant_buckets
    ),
    bucket_mechanisms = bucket_mechanisms,
    summary = list(
      mechanism_split = mechanism_summary$mechanism_split,
      all_representative_abs_att_ratios_below_one = isTRUE(
        mechanism_summary$all_representative_abs_att_ratios_below_one
      ),
      all_buckets_ranked_by_non_near_threshold_mass = isTRUE(
        mechanism_summary$all_buckets_ranked_by_non_near_threshold_mass
      )
    )
  )
}


.story_local_wcb_case_combo_mechanism_bridge <- function() {
  case_frontier <- .story_local_wcb_non_near_threshold_case_frontier()
  mechanism_summary <- .story_local_wcb_representative_far_miss_mechanism_summary()
  minimum_case_rank <- as.integer(case_frontier$minimum_case_rank_to_clear_threshold)
  bucket_frequency <- case_frontier$summary$minimum_case_combo_bucket_frequency
  top_bucket_order <- as.integer(mechanism_summary$top_bucket_order)
  low_signal_plus_elevated_se_buckets <- as.integer(
    mechanism_summary$low_signal_plus_elevated_se_buckets
  )
  low_signal_dominant_buckets <- as.integer(
    mechanism_summary$low_signal_dominant_buckets
  )
  residual_bucket_support <- setdiff(
    as.integer(case_frontier$summary$minimum_case_combo_bucket_support),
    top_bucket_order
  )
  sum_bucket_frequency <- function(bucket_ids) {
    if (length(bucket_ids) == 0L || length(bucket_frequency) == 0L) {
      return(0L)
    }
    values <- bucket_frequency[match(as.character(bucket_ids), names(bucket_frequency))]
    values[is.na(values)] <- 0L
    as.integer(sum(values))
  }

  top_bucket_case_count <- sum_bucket_frequency(top_bucket_order)
  low_signal_plus_case_count <- sum_bucket_frequency(low_signal_plus_elevated_se_buckets)
  low_signal_dominant_case_count <- sum_bucket_frequency(low_signal_dominant_buckets)
  residual_case_count <- sum_bucket_frequency(residual_bucket_support)
  dominant_bucket_case_count <- if (length(bucket_frequency) == 0L) {
    0L
  } else {
    as.integer(unname(bucket_frequency[[1L]]))
  }

  list(
    case_id = case_frontier$case_id,
    exact_status = "story-local-case-combo-mechanism-bridge-frozen",
    numeric_status = "case-combo-mechanism-family-shares-machine-readable",
    blocker_boundary = case_frontier$blocker_boundary,
    story_live_blockers = case_frontier$story_live_blockers,
    consumer_summary_source = paste0(
      "story_local_wcb_non_near_threshold_case_frontier()+",
      "story_local_wcb_representative_far_miss_mechanism_summary()"
    ),
    minimum_case_rank_to_clear_threshold = minimum_case_rank,
    minimum_case_combo_bucket_frequency = bucket_frequency,
    top_bucket_order = top_bucket_order,
    low_signal_plus_elevated_se_buckets = low_signal_plus_elevated_se_buckets,
    low_signal_dominant_buckets = low_signal_dominant_buckets,
    residual_bucket_support = as.integer(residual_bucket_support),
    summary_consistency_checks = list(
      top_bucket_and_residual_reconstruct_case_combo = (
        top_bucket_case_count + residual_case_count
      ) == minimum_case_rank,
      mechanism_split_reconstructs_top_bucket_case_count = (
        low_signal_plus_case_count + low_signal_dominant_case_count
      ) == top_bucket_case_count
    ),
    summary = list(
      top_bucket_case_count = as.integer(top_bucket_case_count),
      top_bucket_case_share = if (minimum_case_rank > 0L) {
        top_bucket_case_count / minimum_case_rank
      } else {
        NaN
      },
      low_signal_plus_elevated_se_case_count = as.integer(low_signal_plus_case_count),
      low_signal_plus_elevated_se_case_share = if (minimum_case_rank > 0L) {
        low_signal_plus_case_count / minimum_case_rank
      } else {
        NaN
      },
      low_signal_dominant_case_count = as.integer(low_signal_dominant_case_count),
      low_signal_dominant_case_share = if (minimum_case_rank > 0L) {
        low_signal_dominant_case_count / minimum_case_rank
      } else {
        NaN
      },
      residual_case_count = as.integer(residual_case_count),
      residual_case_share = if (minimum_case_rank > 0L) {
        residual_case_count / minimum_case_rank
      } else {
        NaN
      },
      dominant_bucket_case_count = as.integer(dominant_bucket_case_count),
      dominant_bucket_share_of_minimum_case_combo = if (minimum_case_rank > 0L) {
        dominant_bucket_case_count / minimum_case_rank
      } else {
        NaN
      },
      top_bucket_family_covers_threshold_combo_majority = (
        top_bucket_case_count * 2L
      ) > minimum_case_rank
    )
  )
}



.story_local_wcb_partition_exclusive_bucket_representatives <- function() {
  partition <- .story_local_wcb_case_combo_residual_partition()
  case_frontier <- .story_local_wcb_non_near_threshold_case_frontier()
  residual_replays <- .story_local_wcb_non_near_threshold_case_frontier_residual_replays()

  threshold_only_bucket_support <- as.integer(partition$minimum_only_bucket_support)
  residual_only_bucket_support <- as.integer(partition$residual_only_bucket_support)
  threshold_pool <- if (length(case_frontier$case_repair_frontier) == 0L) {
    list()
  } else {
    case_frontier$case_repair_frontier[seq_len(partition$minimum_case_rank_to_clear_threshold)]
  }
  threshold_only_representatives <- Filter(
    function(row) row$treated_cluster_count %in% threshold_only_bucket_support,
    threshold_pool
  )
  threshold_only_replay_seeds <- if (length(threshold_only_representatives) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      threshold_only_representatives,
      `[[`,
      integer(1L),
      "replay_seed"
    ))
  }
  threshold_targeted_replays <- if (length(threshold_only_replay_seeds) == 0L) {
    list()
  } else {
    .story_local_wcb_power_diagnostics(
      replay_seeds = threshold_only_replay_seeds
    )$result$targeted_replays
  }
  threshold_only_representatives <- lapply(
    seq_along(threshold_only_representatives),
    function(idx) {
      c(
        threshold_only_representatives[[idx]],
        list(targeted_replay = threshold_targeted_replays[[idx]])
      )
    }
  )

  residual_only_representatives <- Filter(
    function(row) row$treated_cluster_count %in% residual_only_bucket_support,
    residual_replays$residual_case_frontier
  )
  residual_only_replay_seeds <- if (length(residual_only_representatives) == 0L) {
    integer(0L)
  } else {
    as.integer(vapply(
      residual_only_representatives,
      `[[`,
      integer(1L),
      "replay_seed"
    ))
  }

  threshold_only_gap_above_threshold <- if (length(threshold_only_representatives) == 0L) {
    numeric(0L)
  } else {
    vapply(
      threshold_only_representatives,
      `[[`,
      numeric(1L),
      "pvalue_gap_above_threshold"
    )
  }
  threshold_only_t_stats <- if (length(threshold_only_representatives) == 0L) {
    numeric(0L)
  } else {
    vapply(
      lapply(threshold_only_representatives, `[[`, "targeted_replay"),
      `[[`,
      numeric(1L),
      "t_stat_original"
    )
  }
  residual_only_gap_above_threshold <- if (length(residual_only_representatives) == 0L) {
    numeric(0L)
  } else {
    vapply(
      residual_only_representatives,
      `[[`,
      numeric(1L),
      "pvalue_gap_above_threshold"
    )
  }

  closest_threshold_only_case <- if (length(threshold_only_representatives) == 0L) {
    NULL
  } else {
    threshold_only_representatives[[which.min(threshold_only_gap_above_threshold)]]
  }
  highest_t_threshold_only_case <- if (length(threshold_only_representatives) == 0L) {
    NULL
  } else {
    threshold_only_representatives[[which.max(abs(threshold_only_t_stats))]]
  }
  easiest_residual_only_case <- if (length(residual_only_representatives) == 0L) {
    NULL
  } else {
    residual_only_representatives[[which.min(residual_only_gap_above_threshold)]]
  }
  hardest_residual_only_case <- if (length(residual_only_representatives) == 0L) {
    NULL
  } else {
    residual_only_representatives[[which.max(residual_only_gap_above_threshold)]]
  }

  total_partition_exclusive_cases <- length(threshold_only_representatives) +
    length(residual_only_representatives)
  all_partition_exclusive_targeted_replays <- c(
    lapply(threshold_only_representatives, `[[`, "targeted_replay"),
    lapply(residual_only_representatives, `[[`, "targeted_replay")
  )
  residual_only_bucket_reconstruction <- if (length(residual_only_representatives) == 0L) {
    FALSE
  } else {
    identical(
      as.integer(unique(vapply(
        residual_only_representatives,
        `[[`,
        integer(1L),
        "treated_cluster_count"
      ))),
      7L
    )
  }

  list(
    case_id = partition$case_id,
    exact_status =
      "story-local-partition-exclusive-bucket-representatives-frozen",
    numeric_status =
      "threshold-only-vs-residual-only-representative-replays-machine-readable",
    blocker_boundary = partition$blocker_boundary,
    story_live_blockers = partition$story_live_blockers,
    consumer_summary_source = paste0(
      "story_local_wcb_case_combo_residual_partition()+",
      "story_local_wcb_power_diagnostics(replay_seeds = partition_exclusive_bucket_representatives)"
    ),
    threshold_only_bucket_support = threshold_only_bucket_support,
    residual_only_bucket_support = if (length(residual_only_bucket_support) == 1L) {
      residual_only_bucket_support[[1L]]
    } else {
      residual_only_bucket_support
    },
    threshold_only_replay_seeds = threshold_only_replay_seeds,
    residual_only_replay_seeds = residual_only_replay_seeds,
    threshold_only_representatives = threshold_only_representatives,
    residual_only_representatives = residual_only_representatives,
    summary = list(
      threshold_only_case_count = as.integer(length(threshold_only_representatives)),
      residual_only_case_count = as.integer(length(residual_only_representatives)),
      threshold_only_share_of_partition_exclusive_cases = if (
        total_partition_exclusive_cases > 0L
      ) {
        length(threshold_only_representatives) / total_partition_exclusive_cases
      } else {
        NaN
      },
      residual_only_share_of_partition_exclusive_cases = if (
        total_partition_exclusive_cases > 0L
      ) {
        length(residual_only_representatives) / total_partition_exclusive_cases
      } else {
        NaN
      },
      closest_threshold_only_replay_seed = if (is.null(closest_threshold_only_case)) {
        NULL
      } else {
        as.integer(closest_threshold_only_case$replay_seed)
      },
      closest_threshold_only_gap_above_threshold = if (is.null(closest_threshold_only_case)) {
        NULL
      } else {
        as.numeric(closest_threshold_only_case$pvalue_gap_above_threshold)
      },
      highest_threshold_only_t_stat_seed = if (is.null(highest_t_threshold_only_case)) {
        NULL
      } else {
        as.integer(highest_t_threshold_only_case$replay_seed)
      },
      highest_threshold_only_t_stat = if (is.null(highest_t_threshold_only_case)) {
        NULL
      } else {
        as.numeric(abs(highest_t_threshold_only_case$targeted_replay$t_stat_original))
      },
      easiest_residual_only_replay_seed = if (is.null(easiest_residual_only_case)) {
        NULL
      } else {
        as.integer(easiest_residual_only_case$replay_seed)
      },
      hardest_residual_only_replay_seed = if (is.null(hardest_residual_only_case)) {
        NULL
      } else {
        as.integer(hardest_residual_only_case$replay_seed)
      },
      all_partition_exclusive_targeted_replays_remain_misses = all(!vapply(
        all_partition_exclusive_targeted_replays,
        `[[`,
        logical(1L),
        "reject"
      )),
      bucket8_is_closest_threshold_only_case = !is.null(closest_threshold_only_case) &&
        identical(as.integer(closest_threshold_only_case$treated_cluster_count), 8L),
      bucket9_has_highest_threshold_only_t_stat =
        !is.null(highest_t_threshold_only_case) &&
          identical(
            as.integer(highest_t_threshold_only_case$treated_cluster_count),
            9L
          ),
      bucket7_reconstructs_residual_only_family = residual_only_bucket_reconstruction
    )
  )
}


.story_local_wcb_partition_exclusive_bucket_signal_contrast <- function() {
  exclusive <- .story_local_wcb_partition_exclusive_bucket_representatives()
  threshold_only_representatives <- exclusive$threshold_only_representatives
  residual_only_representatives <- exclusive$residual_only_representatives

  select_closest_threshold_only_case <- if (length(threshold_only_representatives) == 0L) {
    NULL
  } else {
    threshold_only_representatives[[which.min(vapply(
      threshold_only_representatives,
      `[[`,
      numeric(1L),
      "pvalue_gap_above_threshold"
    ))]]
  }
  select_highest_t_threshold_only_case <- if (length(threshold_only_representatives) == 0L) {
    NULL
  } else {
    threshold_only_representatives[[which.max(vapply(
      lapply(threshold_only_representatives, `[[`, "targeted_replay"),
      function(row) abs(as.numeric(row$t_stat_original)),
      numeric(1L)
    ))]]
  }
  select_easiest_residual_only_case <- if (length(residual_only_representatives) == 0L) {
    NULL
  } else {
    residual_only_representatives[[which.min(vapply(
      residual_only_representatives,
      `[[`,
      numeric(1L),
      "pvalue_gap_above_threshold"
    ))]]
  }
  select_hardest_residual_only_case <- if (length(residual_only_representatives) == 0L) {
    NULL
  } else {
    residual_only_representatives[[which.max(vapply(
      residual_only_representatives,
      `[[`,
      numeric(1L),
      "pvalue_gap_above_threshold"
    ))]]
  }

  compact_row <- function(row) {
    if (is.null(row)) {
      return(NULL)
    }

    list(
      case_rank = as.integer(row$case_rank),
      replay_seed = as.integer(row$replay_seed),
      replay_seed_offset = as.integer(row$replay_seed_offset),
      replay_attempt_id = as.integer(row$attempt_id),
      treated_cluster_count = as.integer(row$treated_cluster_count),
      pvalue_gap_above_threshold = as.numeric(row$pvalue_gap_above_threshold),
      targeted_replay = list(
        seed = as.integer(row$targeted_replay$seed),
        reject = isTRUE(row$targeted_replay$reject),
        pvalue = as.numeric(row$targeted_replay$pvalue),
        att = as.numeric(row$targeted_replay$att),
        original_se = as.numeric(row$targeted_replay$original_se),
        t_stat_original = as.numeric(row$targeted_replay$t_stat_original)
      )
    )
  }
  abs_att <- function(row) {
    if (is.null(row)) {
      return(NA_real_)
    }
    abs(as.numeric(row$targeted_replay$att))
  }
  original_se <- function(row) {
    if (is.null(row)) {
      return(NA_real_)
    }
    as.numeric(row$targeted_replay$original_se)
  }
  abs_t_stat <- function(row) {
    if (is.null(row)) {
      return(NA_real_)
    }
    abs(as.numeric(row$targeted_replay$t_stat_original))
  }
  pvalue_gap <- function(row) {
    if (is.null(row)) {
      return(NA_real_)
    }
    as.numeric(row$pvalue_gap_above_threshold)
  }

  threshold_abs_att <- if (length(threshold_only_representatives) == 0L) {
    numeric(0L)
  } else {
    vapply(threshold_only_representatives, abs_att, numeric(1L))
  }
  residual_abs_att <- if (length(residual_only_representatives) == 0L) {
    numeric(0L)
  } else {
    vapply(residual_only_representatives, abs_att, numeric(1L))
  }
  threshold_abs_t <- if (length(threshold_only_representatives) == 0L) {
    numeric(0L)
  } else {
    vapply(threshold_only_representatives, abs_t_stat, numeric(1L))
  }
  residual_abs_t <- if (length(residual_only_representatives) == 0L) {
    numeric(0L)
  } else {
    vapply(residual_only_representatives, abs_t_stat, numeric(1L))
  }

  bucket9_abs_att_gain_over_bucket8 <- abs_att(select_highest_t_threshold_only_case) -
    abs_att(select_closest_threshold_only_case)
  bucket9_original_se_reduction_vs_bucket8 <- original_se(select_closest_threshold_only_case) -
    original_se(select_highest_t_threshold_only_case)
  bucket9_t_stat_gain_over_bucket8 <- abs_t_stat(select_highest_t_threshold_only_case) -
    abs_t_stat(select_closest_threshold_only_case)
  bucket9_pvalue_gap_penalty_vs_bucket8 <- pvalue_gap(select_highest_t_threshold_only_case) -
    pvalue_gap(select_closest_threshold_only_case)
  bucket8_abs_att_gain_over_easiest_residual <- abs_att(select_closest_threshold_only_case) -
    abs_att(select_easiest_residual_only_case)
  bucket8_original_se_penalty_vs_easiest_residual <- original_se(select_closest_threshold_only_case) -
    original_se(select_easiest_residual_only_case)
  bucket8_t_stat_gain_over_easiest_residual <- abs_t_stat(select_closest_threshold_only_case) -
    abs_t_stat(select_easiest_residual_only_case)
  bucket9_abs_att_gain_over_hardest_residual <- abs_att(select_highest_t_threshold_only_case) -
    abs_att(select_hardest_residual_only_case)
  bucket9_original_se_advantage_over_hardest_residual <- original_se(select_hardest_residual_only_case) -
    original_se(select_highest_t_threshold_only_case)
  bucket9_t_stat_gain_over_hardest_residual <- abs_t_stat(select_highest_t_threshold_only_case) -
    abs_t_stat(select_hardest_residual_only_case)

  list(
    case_id = exclusive$case_id,
    exact_status = "story-local-partition-exclusive-bucket-signal-contrast-frozen",
    numeric_status = "threshold-only-vs-residual-only-signal-se-contrast-machine-readable",
    blocker_boundary = exclusive$blocker_boundary,
    story_live_blockers = exclusive$story_live_blockers,
    consumer_summary_source = "story_local_wcb_partition_exclusive_bucket_representatives()",
    threshold_only_replay_seeds = as.integer(exclusive$threshold_only_replay_seeds),
    residual_only_replay_seeds = as.integer(exclusive$residual_only_replay_seeds),
    threshold_only_pair = list(
      closest_threshold_only = compact_row(select_closest_threshold_only_case),
      highest_t_threshold_only = compact_row(select_highest_t_threshold_only_case)
    ),
    residual_only_pair = list(
      easiest_residual_only = compact_row(select_easiest_residual_only_case),
      hardest_residual_only = compact_row(select_hardest_residual_only_case)
    ),
    summary = list(
      bucket9_abs_att_gain_over_bucket8 = as.numeric(bucket9_abs_att_gain_over_bucket8),
      bucket9_original_se_reduction_vs_bucket8 = as.numeric(
        bucket9_original_se_reduction_vs_bucket8
      ),
      bucket9_t_stat_gain_over_bucket8 = as.numeric(bucket9_t_stat_gain_over_bucket8),
      bucket9_pvalue_gap_penalty_vs_bucket8 = as.numeric(
        bucket9_pvalue_gap_penalty_vs_bucket8
      ),
      bucket8_abs_att_gain_over_easiest_residual = as.numeric(
        bucket8_abs_att_gain_over_easiest_residual
      ),
      bucket8_original_se_penalty_vs_easiest_residual = as.numeric(
        bucket8_original_se_penalty_vs_easiest_residual
      ),
      bucket8_t_stat_gain_over_easiest_residual = as.numeric(
        bucket8_t_stat_gain_over_easiest_residual
      ),
      bucket9_abs_att_gain_over_hardest_residual = as.numeric(
        bucket9_abs_att_gain_over_hardest_residual
      ),
      bucket9_original_se_advantage_over_hardest_residual = as.numeric(
        bucket9_original_se_advantage_over_hardest_residual
      ),
      bucket9_t_stat_gain_over_hardest_residual = as.numeric(
        bucket9_t_stat_gain_over_hardest_residual
      ),
      all_threshold_only_abs_att_exceed_residual_only_family = (
        length(threshold_abs_att) > 0L &&
          length(residual_abs_att) > 0L &&
          all(threshold_abs_att > max(residual_abs_att))
      ),
      all_threshold_only_t_stats_exceed_residual_only_family = (
        length(threshold_abs_t) > 0L &&
          length(residual_abs_t) > 0L &&
          all(threshold_abs_t > max(residual_abs_t))
      ),
      bucket9_pairs_more_signal_with_lower_se_than_bucket8 = (
        is.finite(bucket9_abs_att_gain_over_bucket8) &&
          is.finite(bucket9_original_se_reduction_vs_bucket8) &&
          is.finite(bucket9_t_stat_gain_over_bucket8) &&
          bucket9_abs_att_gain_over_bucket8 > 0 &&
          bucket9_original_se_reduction_vs_bucket8 > 0 &&
          bucket9_t_stat_gain_over_bucket8 > 0
      ),
      bucket8_retains_higher_t_than_easiest_residual_despite_se_penalty = (
        is.finite(bucket8_abs_att_gain_over_easiest_residual) &&
          is.finite(bucket8_original_se_penalty_vs_easiest_residual) &&
          is.finite(bucket8_t_stat_gain_over_easiest_residual) &&
          bucket8_abs_att_gain_over_easiest_residual > 0 &&
          bucket8_original_se_penalty_vs_easiest_residual > 0 &&
          bucket8_t_stat_gain_over_easiest_residual > 0
      ),
      bucket9_dominates_hardest_residual_on_signal_and_precision = (
        is.finite(bucket9_abs_att_gain_over_hardest_residual) &&
          is.finite(bucket9_original_se_advantage_over_hardest_residual) &&
          is.finite(bucket9_t_stat_gain_over_hardest_residual) &&
          bucket9_abs_att_gain_over_hardest_residual > 0 &&
          bucket9_original_se_advantage_over_hardest_residual > 0 &&
          bucket9_t_stat_gain_over_hardest_residual > 0
      )
    )
  )
}


.story_local_wcb_partition_exclusive_residual_monotone_decay <- function() {
  contrast <- .story_local_wcb_partition_exclusive_bucket_signal_contrast()
  exclusive <- .story_local_wcb_partition_exclusive_bucket_representatives()
  residual_only_representatives <- exclusive$residual_only_representatives

  ordered_residual_only_cases <- if (length(residual_only_representatives) == 0L) {
    list()
  } else {
    residual_only_representatives[order(vapply(
      residual_only_representatives,
      `[[`,
      numeric(1L),
      "pvalue_gap_above_threshold"
    ))]
  }

  select_easiest_residual_only_case <- if (length(ordered_residual_only_cases) >= 1L) {
    ordered_residual_only_cases[[1L]]
  } else {
    NULL
  }
  select_middle_residual_only_case <- if (length(ordered_residual_only_cases) >= 2L) {
    ordered_residual_only_cases[[2L]]
  } else {
    NULL
  }
  select_hardest_residual_only_case <- if (length(ordered_residual_only_cases) >= 3L) {
    ordered_residual_only_cases[[3L]]
  } else {
    NULL
  }

  compact_row <- function(row) {
    if (is.null(row)) {
      return(NULL)
    }

    list(
      case_rank = as.integer(row$case_rank),
      replay_seed = as.integer(row$replay_seed),
      replay_seed_offset = as.integer(row$replay_seed_offset),
      replay_attempt_id = as.integer(row$attempt_id),
      treated_cluster_count = as.integer(row$treated_cluster_count),
      pvalue_gap_above_threshold = as.numeric(row$pvalue_gap_above_threshold),
      targeted_replay = list(
        seed = as.integer(row$targeted_replay$seed),
        reject = isTRUE(row$targeted_replay$reject),
        pvalue = as.numeric(row$targeted_replay$pvalue),
        att = as.numeric(row$targeted_replay$att),
        original_se = as.numeric(row$targeted_replay$original_se),
        t_stat_original = as.numeric(row$targeted_replay$t_stat_original)
      )
    )
  }
  abs_att <- function(row) {
    if (is.null(row)) {
      return(NA_real_)
    }
    abs(as.numeric(row$targeted_replay$att))
  }
  original_se <- function(row) {
    if (is.null(row)) {
      return(NA_real_)
    }
    as.numeric(row$targeted_replay$original_se)
  }
  abs_t_stat <- function(row) {
    if (is.null(row)) {
      return(NA_real_)
    }
    abs(as.numeric(row$targeted_replay$t_stat_original))
  }
  pvalue_gap <- function(row) {
    if (is.null(row)) {
      return(NA_real_)
    }
    as.numeric(row$pvalue_gap_above_threshold)
  }
  safe_ratio <- function(numerator, denominator) {
    if (!is.finite(numerator) || !is.finite(denominator) || denominator == 0) {
      return(NA_real_)
    }
    as.numeric(numerator / denominator)
  }

  bucket48_abs_att_gain_over_bucket53 <- abs_att(select_easiest_residual_only_case) -
    abs_att(select_middle_residual_only_case)
  bucket53_abs_att_gain_over_bucket79 <- abs_att(select_middle_residual_only_case) -
    abs_att(select_hardest_residual_only_case)
  bucket53_abs_att_retention_vs_bucket48 <- safe_ratio(
    abs_att(select_middle_residual_only_case),
    abs_att(select_easiest_residual_only_case)
  )
  bucket79_abs_att_retention_vs_bucket53 <- safe_ratio(
    abs_att(select_hardest_residual_only_case),
    abs_att(select_middle_residual_only_case)
  )
  bucket53_original_se_penalty_vs_bucket48 <- original_se(select_middle_residual_only_case) -
    original_se(select_easiest_residual_only_case)
  bucket79_original_se_penalty_vs_bucket53 <- original_se(select_hardest_residual_only_case) -
    original_se(select_middle_residual_only_case)
  bucket53_original_se_multiplier_vs_bucket48 <- safe_ratio(
    original_se(select_middle_residual_only_case),
    original_se(select_easiest_residual_only_case)
  )
  bucket79_original_se_multiplier_vs_bucket53 <- safe_ratio(
    original_se(select_hardest_residual_only_case),
    original_se(select_middle_residual_only_case)
  )
  bucket48_t_stat_gain_over_bucket53 <- abs_t_stat(select_easiest_residual_only_case) -
    abs_t_stat(select_middle_residual_only_case)
  bucket53_t_stat_gain_over_bucket79 <- abs_t_stat(select_middle_residual_only_case) -
    abs_t_stat(select_hardest_residual_only_case)
  bucket53_t_stat_retention_vs_bucket48 <- safe_ratio(
    abs_t_stat(select_middle_residual_only_case),
    abs_t_stat(select_easiest_residual_only_case)
  )
  bucket79_t_stat_retention_vs_bucket53 <- safe_ratio(
    abs_t_stat(select_hardest_residual_only_case),
    abs_t_stat(select_middle_residual_only_case)
  )
  bucket53_pvalue_gap_penalty_vs_bucket48 <- pvalue_gap(select_middle_residual_only_case) -
    pvalue_gap(select_easiest_residual_only_case)
  bucket79_pvalue_gap_penalty_vs_bucket53 <- pvalue_gap(select_hardest_residual_only_case) -
    pvalue_gap(select_middle_residual_only_case)
  bucket53_pvalue_gap_multiplier_vs_bucket48 <- safe_ratio(
    pvalue_gap(select_middle_residual_only_case),
    pvalue_gap(select_easiest_residual_only_case)
  )
  bucket79_pvalue_gap_multiplier_vs_bucket53 <- safe_ratio(
    pvalue_gap(select_hardest_residual_only_case),
    pvalue_gap(select_middle_residual_only_case)
  )

  residual_abs_att <- vapply(ordered_residual_only_cases, abs_att, numeric(1L))
  residual_abs_t <- vapply(ordered_residual_only_cases, abs_t_stat, numeric(1L))
  residual_original_se <- vapply(ordered_residual_only_cases, original_se, numeric(1L))
  residual_pvalue_gap <- vapply(ordered_residual_only_cases, pvalue_gap, numeric(1L))
  residual_case_ranks <- vapply(
    ordered_residual_only_cases,
    `[[`,
    integer(1L),
    "case_rank"
  )
  residual_replay_seed_offsets <- vapply(
    ordered_residual_only_cases,
    `[[`,
    integer(1L),
    "replay_seed_offset"
  )
  residual_replay_attempt_ids <- vapply(
    ordered_residual_only_cases,
    `[[`,
    integer(1L),
    "attempt_id"
  )

  list(
    case_id = contrast$case_id,
    exact_status = "story-local-partition-exclusive-residual-monotone-decay-frozen",
    numeric_status = "residual-only-48-53-79-monotone-decay-machine-readable",
    blocker_boundary = contrast$blocker_boundary,
    story_live_blockers = contrast$story_live_blockers,
    consumer_summary_source = "story_local_wcb_partition_exclusive_bucket_signal_contrast()",
    threshold_only_replay_seeds = as.integer(contrast$threshold_only_replay_seeds),
    residual_only_replay_seeds = as.integer(exclusive$residual_only_replay_seeds),
    threshold_only_anchor_pair = list(
      closest_threshold_only = compact_row(
        contrast$threshold_only_pair$closest_threshold_only
      ),
      highest_t_threshold_only = compact_row(
        contrast$threshold_only_pair$highest_t_threshold_only
      )
    ),
    residual_only_triplet = list(
      easiest_residual_only = compact_row(select_easiest_residual_only_case),
      middle_residual_only = compact_row(select_middle_residual_only_case),
      hardest_residual_only = compact_row(select_hardest_residual_only_case)
    ),
    summary = list(
      bucket48_abs_att_gain_over_bucket53 = as.numeric(
        bucket48_abs_att_gain_over_bucket53
      ),
      bucket53_abs_att_gain_over_bucket79 = as.numeric(
        bucket53_abs_att_gain_over_bucket79
      ),
      bucket53_abs_att_retention_vs_bucket48 = as.numeric(
        bucket53_abs_att_retention_vs_bucket48
      ),
      bucket79_abs_att_retention_vs_bucket53 = as.numeric(
        bucket79_abs_att_retention_vs_bucket53
      ),
      bucket53_original_se_penalty_vs_bucket48 = as.numeric(
        bucket53_original_se_penalty_vs_bucket48
      ),
      bucket79_original_se_penalty_vs_bucket53 = as.numeric(
        bucket79_original_se_penalty_vs_bucket53
      ),
      bucket53_original_se_multiplier_vs_bucket48 = as.numeric(
        bucket53_original_se_multiplier_vs_bucket48
      ),
      bucket79_original_se_multiplier_vs_bucket53 = as.numeric(
        bucket79_original_se_multiplier_vs_bucket53
      ),
      bucket48_t_stat_gain_over_bucket53 = as.numeric(
        bucket48_t_stat_gain_over_bucket53
      ),
      bucket53_t_stat_gain_over_bucket79 = as.numeric(
        bucket53_t_stat_gain_over_bucket79
      ),
      bucket53_t_stat_retention_vs_bucket48 = as.numeric(
        bucket53_t_stat_retention_vs_bucket48
      ),
      bucket79_t_stat_retention_vs_bucket53 = as.numeric(
        bucket79_t_stat_retention_vs_bucket53
      ),
      bucket53_pvalue_gap_penalty_vs_bucket48 = as.numeric(
        bucket53_pvalue_gap_penalty_vs_bucket48
      ),
      bucket79_pvalue_gap_penalty_vs_bucket53 = as.numeric(
        bucket79_pvalue_gap_penalty_vs_bucket53
      ),
      bucket53_pvalue_gap_multiplier_vs_bucket48 = as.numeric(
        bucket53_pvalue_gap_multiplier_vs_bucket48
      ),
      bucket79_pvalue_gap_multiplier_vs_bucket53 = as.numeric(
        bucket79_pvalue_gap_multiplier_vs_bucket53
      ),
      residual_only_abs_att_strictly_descends = (
        length(residual_abs_att) == 3L &&
          all(diff(residual_abs_att) < 0)
      ),
      residual_only_t_stats_strictly_descend = (
        length(residual_abs_t) == 3L &&
          all(diff(residual_abs_t) < 0)
      ),
      residual_only_original_se_strictly_ascends = (
        length(residual_original_se) == 3L &&
          all(diff(residual_original_se) > 0)
      ),
      residual_only_gap_above_threshold_strictly_ascends = (
        length(residual_pvalue_gap) == 3L &&
          all(diff(residual_pvalue_gap) > 0)
      ),
      case_ranks_strictly_ascend = (
        length(residual_case_ranks) == 3L &&
          all(diff(residual_case_ranks) > 0L)
      ),
      replay_seed_offsets_strictly_ascend = (
        length(residual_replay_seed_offsets) == 3L &&
          all(diff(residual_replay_seed_offsets) > 0L)
      ),
      replay_attempt_ids_strictly_ascend = (
        length(residual_replay_attempt_ids) == 3L &&
          all(diff(residual_replay_attempt_ids) > 0L)
      ),
      threshold_only_anchor_abs_att_still_exceeds_residual_family = isTRUE(
        contrast$summary$all_threshold_only_abs_att_exceed_residual_only_family
      ),
      threshold_only_anchor_t_stats_still_exceed_residual_family = isTRUE(
        contrast$summary$all_threshold_only_t_stats_exceed_residual_only_family
      )
    )
  )
}


.story_local_wcb_partition_exclusive_residual_targeted_replay_rank_stability <- function() {
  exclusive_decay <- .story_local_wcb_partition_exclusive_residual_monotone_decay()
  replay_diagnostics <- .story_local_wcb_power_diagnostics(
    replay_seeds = c(48L, 53L, 79L)
  )
  replay_triplet <- replay_diagnostics$result$targeted_replays
  replay_summary <- replay_diagnostics$result$targeted_replay_summary
  helper_triplet <- exclusive_decay$residual_only_triplet

  compact_replay_triplet_row <- function(helper_row, replay_row) {
    list(
      case_rank = as.integer(helper_row$case_rank),
      replay_seed = as.integer(helper_row$replay_seed),
      replay_seed_offset = as.integer(replay_row$replay_seed_offset),
      replay_attempt_id = as.integer(replay_row$attempt_id),
      replication_id = as.integer(replay_row$replication_id),
      treated_cluster_count = as.integer(replay_row$treated_cluster_count),
      reject = isTRUE(replay_row$reject),
      pvalue = as.numeric(replay_row$pvalue),
      abs_gap_to_threshold = as.numeric(replay_row$abs_gap_to_threshold),
      att = as.numeric(replay_row$att),
      original_se = as.numeric(replay_row$original_se),
      t_stat_original = as.numeric(replay_row$t_stat_original)
    )
  }

  triplet_rows <- lapply(
    c("easiest_residual_only", "middle_residual_only", "hardest_residual_only"),
    function(name) {
      helper_row <- helper_triplet[[name]]
      matching_index <- which(vapply(
        replay_triplet,
        function(row) {
          identical(as.integer(row$seed), as.integer(helper_row$replay_seed))
        },
        logical(1L)
      ))
      if (!length(matching_index)) {
        stop(
          sprintf("Replay diagnostics missing seed %s", helper_row$replay_seed),
          call. = FALSE
        )
      }
      compact_replay_triplet_row(helper_row, replay_triplet[[matching_index[[1L]]]])
    }
  )

  abs_att_values <- unname(as.numeric(vapply(
    triplet_rows,
    function(row) abs(row$att),
    numeric(1L)
  )))
  original_se_values <- unname(as.numeric(vapply(
    triplet_rows,
    `[[`,
    numeric(1L),
    "original_se"
  )))
  t_stat_values <- unname(as.numeric(vapply(
    triplet_rows,
    function(row) abs(row$t_stat_original),
    numeric(1L)
  )))
  pvalue_values <- unname(as.numeric(vapply(
    triplet_rows,
    `[[`,
    numeric(1L),
    "pvalue"
  )))
  pvalue_gap_values <- unname(as.numeric(vapply(
    triplet_rows,
    `[[`,
    numeric(1L),
    "abs_gap_to_threshold"
  )))
  replay_seed_offsets <- unname(as.integer(vapply(
    triplet_rows,
    `[[`,
    integer(1L),
    "replay_seed_offset"
  )))
  replay_attempt_ids <- unname(as.integer(vapply(
    triplet_rows,
    `[[`,
    integer(1L),
    "replay_attempt_id"
  )))
  replay_case_ranks <- unname(as.integer(vapply(
    triplet_rows,
    `[[`,
    integer(1L),
    "case_rank"
  )))
  replication_ids <- unname(as.integer(vapply(
    triplet_rows,
    `[[`,
    integer(1L),
    "replication_id"
  )))
  treated_cluster_pattern <- unname(as.integer(vapply(
    triplet_rows,
    `[[`,
    integer(1L),
    "treated_cluster_count"
  )))
  helper_rejection_pattern <- unname(as.logical(vapply(
    triplet_rows,
    `[[`,
    logical(1L),
    "reject"
  )))
  triplet_tuples <- lapply(
    triplet_rows,
    function(row) {
      list(
        replay_seed = as.integer(row$replay_seed),
        replay_attempt_id = as.integer(row$replay_attempt_id),
        replay_case_rank = as.integer(row$case_rank),
        replay_seed_offset = as.integer(row$replay_seed_offset),
        replication_id = as.integer(row$replication_id),
        treated_cluster_count = as.integer(row$treated_cluster_count),
        reject = isTRUE(row$reject),
        pvalue = as.numeric(row$pvalue),
        abs_gap_to_threshold = as.numeric(row$abs_gap_to_threshold),
        abs_att = abs(as.numeric(row$att)),
        original_se = as.numeric(row$original_se),
        t_stat = abs(as.numeric(row$t_stat_original))
      )
    }
  )
  triplet_matches_residual_monotone_helper <- isTRUE(all.equal(
    unname(as.integer(vapply(triplet_rows, `[[`, integer(1L), "replay_seed"))),
    unname(as.integer(exclusive_decay$residual_only_replay_seeds))
  )) && isTRUE(all.equal(
    abs_att_values,
    unname(as.numeric(c(
      abs(helper_triplet$easiest_residual_only$targeted_replay$att),
      abs(helper_triplet$middle_residual_only$targeted_replay$att),
      abs(helper_triplet$hardest_residual_only$targeted_replay$att)
    ))),
    tolerance = 1e-12
  )) && isTRUE(all.equal(
    original_se_values,
    unname(as.numeric(c(
      helper_triplet$easiest_residual_only$targeted_replay$original_se,
      helper_triplet$middle_residual_only$targeted_replay$original_se,
      helper_triplet$hardest_residual_only$targeted_replay$original_se
    ))),
    tolerance = 1e-12
  )) && isTRUE(all.equal(
    t_stat_values,
    unname(as.numeric(c(
      abs(helper_triplet$easiest_residual_only$targeted_replay$t_stat_original),
      abs(helper_triplet$middle_residual_only$targeted_replay$t_stat_original),
      abs(helper_triplet$hardest_residual_only$targeted_replay$t_stat_original)
    ))),
    tolerance = 1e-12
  ))

  list(
    case_id = exclusive_decay$case_id,
    exact_status = "story-local-partition-exclusive-residual-targeted-replay-rank-stability-frozen",
    numeric_status = "residual-only-48-53-79-targeted-replay-rank-stability-machine-readable",
    blocker_boundary = exclusive_decay$blocker_boundary,
    story_live_blockers = exclusive_decay$story_live_blockers,
    consumer_summary_source = paste0(
      "story_local_wcb_partition_exclusive_residual_monotone_decay()+",
      "story_local_wcb_power_diagnostics(replay_seeds = c(48, 53, 79))"
    ),
    residual_only_replay_seeds = unname(as.integer(vapply(
      triplet_rows,
      `[[`,
      integer(1L),
      "replay_seed"
    ))),
    replay_seed_offsets = replay_seed_offsets,
    replay_attempt_ids = replay_attempt_ids,
    replay_case_ranks = replay_case_ranks,
    treated_cluster_pattern = treated_cluster_pattern,
    helper_rejection_pattern = helper_rejection_pattern,
    weakest_miss = list(
      seed = as.integer(replay_summary$weakest_miss$seed),
      attempt_id = as.integer(replay_summary$weakest_miss$attempt_id),
      replication_id = as.integer(replay_summary$weakest_miss$replication_id),
      replay_seed_offset = as.integer(replay_summary$weakest_miss$replay_seed_offset),
      treated_cluster_count = as.integer(replay_summary$weakest_miss$treated_cluster_count),
      reject = isTRUE(replay_summary$weakest_miss$reject),
      pvalue = as.numeric(replay_summary$weakest_miss$pvalue),
      abs_gap_to_threshold = as.numeric(replay_summary$weakest_miss$abs_gap_to_threshold)
    ),
    triplet_rows = triplet_rows,
    triplet_tuples = triplet_tuples,
    summary = list(
      targeted_replay_count = as.integer(replay_summary$targeted_replay_count),
      replay_seed_labels = as.integer(replay_summary$replay_seed_labels),
      replication_ids = replication_ids,
      abs_att_values = abs_att_values,
      original_se_values = original_se_values,
      t_stat_values = t_stat_values,
      pvalue_values = pvalue_values,
      pvalue_gap_values = pvalue_gap_values,
      replay_attempt_ids_strictly_ascend = (
        length(replay_attempt_ids) == 3L &&
          all(diff(replay_attempt_ids) > 0L)
      ),
      replay_seed_offsets_strictly_ascend = (
        length(replay_seed_offsets) == 3L &&
          all(diff(replay_seed_offsets) > 0L)
      ),
      case_ranks_strictly_ascend = (
        length(replay_case_ranks) == 3L &&
          all(diff(replay_case_ranks) > 0L)
      ),
      abs_att_strictly_descends = (
        length(abs_att_values) == 3L &&
          all(diff(abs_att_values) < 0)
      ),
      original_se_strictly_ascends = (
        length(original_se_values) == 3L &&
          all(diff(original_se_values) > 0)
      ),
      t_stat_strictly_descends = (
        length(t_stat_values) == 3L &&
          all(diff(t_stat_values) < 0)
      ),
      pvalue_strictly_ascends = (
        length(pvalue_values) == 3L &&
          all(diff(pvalue_values) > 0)
      ),
      pvalue_gap_strictly_ascends = (
        length(pvalue_gap_values) == 3L &&
          all(diff(pvalue_gap_values) > 0)
      ),
      all_targeted_replays_remain_non_rejections = !any(helper_rejection_pattern),
      triplet_matches_residual_monotone_helper = triplet_matches_residual_monotone_helper
    )
  )
}


.run_clustering_monte_carlo_relative_gap_scenario <- function(
    scenario, shared_dgp = NULL
) {
  set.seed(as.integer(scenario[["seed"]]))

  n_target <- as.integer(scenario[["n_simulations"]])
  alpha <- if (!is.null(scenario[["alpha"]])) {
    as.numeric(scenario[["alpha"]])
  } else {
    0.05
  }

  wild_spec <- scenario[["estimators"]][["wild"]]
  requested_n_bootstrap <- if (!is.null(wild_spec[["requested_n_bootstrap"]])) {
    as.integer(wild_spec[["requested_n_bootstrap"]])
  } else {
    as.integer(wild_spec[["n_bootstrap"]])
  }

  standard_hits <- logical(n_target)
  wild_hits <- logical(n_target)
  valid_replications <- 0L
  attempts <- 0L
  max_attempts <- n_target * 20L
  last_result <- NULL

  while (valid_replications < n_target) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop(
        sprintf(
          "Clustering Monte Carlo scenario '%s' exceeded %d attempts.",
          scenario[["id"]],
          max_attempts
        ),
        call. = FALSE
      )
    }

    sim_data <- .simulate_clustering_monte_carlo_panel(scenario, shared_dgp)
    standard_result <- .evaluate_clustering_monte_carlo_standard_ci(
      data = sim_data,
      true_tau = as.numeric(scenario[["true_tau"]]),
      alpha = alpha
    )
    if (is.null(standard_result)) {
      next
    }

    wild_result <- .run_wild_cluster_bootstrap_ci(
      data = sim_data,
      requested_n_bootstrap = requested_n_bootstrap,
      weight_type = wild_spec[["weight_type"]],
      alpha = alpha,
      impose_null = wild_spec[["impose_null"]]
    )
    if (!is.finite(wild_result$ci_lower) || !is.finite(wild_result$ci_upper)) {
      next
    }

    valid_replications <- valid_replications + 1L
    standard_hits[[valid_replications]] <- standard_result$covers
    wild_hits[[valid_replications]] <- isTRUE(
      wild_result$ci_lower <= as.numeric(scenario[["true_tau"]]) &&
        as.numeric(scenario[["true_tau"]]) <= wild_result$ci_upper
    )
    last_result <- wild_result
  }

  coverage_std_rate <- mean(standard_hits)
  coverage_wild_rate <- mean(wild_hits)
  coverage_rate_gap <- coverage_wild_rate - coverage_std_rate

  list(
    scenario_id = scenario[["id"]],
    n_simulations = n_target,
    n_attempts = attempts,
    metric_name = "coverage_rate_gap",
    metric_value = coverage_rate_gap,
    metrics = list(
      coverage_std_rate = coverage_std_rate,
      coverage_wild_rate = coverage_wild_rate,
      coverage_rate_gap = coverage_rate_gap
    ),
    requested_n_bootstrap = last_result$requested_n_bootstrap,
    actual_n_bootstrap = last_result$actual_n_bootstrap,
    weight_type = last_result$weight_type,
    impose_null = last_result$impose_null,
    full_enumeration = last_result$full_enumeration,
    ci_method = last_result$ci_method
  )
}


.run_clustering_monte_carlo_scenario <- function(scenario, shared_dgp = NULL) {
  estimator <- scenario[["estimator"]]

  if (identical(estimator, "wild_cluster_bootstrap")) {
    if (!is.null(scenario[["decision_direction"]])) {
      return(
        .run_clustering_monte_carlo_wild_bootstrap_decision_scenario(
          scenario,
          shared_dgp = shared_dgp
        )
      )
    }

    return(
      .run_clustering_monte_carlo_wild_bootstrap_scenario(
        scenario,
        shared_dgp = shared_dgp
      )
    )
  }

  if (!is.null(scenario[["estimators"]])) {
    return(
      .run_clustering_monte_carlo_relative_gap_scenario(
        scenario,
        shared_dgp = shared_dgp
      )
    )
  }

  if (!estimator %in% c("cluster_robust_t_ci", "cluster_robust_t_test")) {
    stop(
      sprintf(
        "Unsupported clustering Monte Carlo estimator '%s' for scenario '%s'.",
        estimator,
        scenario[["id"]]
      ),
      call. = FALSE
    )
  }

  set.seed(as.integer(scenario[["seed"]]))

  n_target <- as.integer(scenario[["n_simulations"]])
  alpha <- if (!is.null(scenario[["alpha"]])) {
    as.numeric(scenario[["alpha"]])
  } else {
    0.05
  }

  metric_hits <- logical(n_target)
  valid_replications <- 0L
  attempts <- 0L
  max_attempts <- n_target * 20L

  while (valid_replications < n_target) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop(
        sprintf(
          "Clustering Monte Carlo scenario '%s' exceeded %d attempts.",
          scenario[["id"]],
          max_attempts
        ),
        call. = FALSE
      )
    }

    sim_data <- .simulate_clustering_monte_carlo_panel(scenario, shared_dgp)
    standard_result <- .evaluate_clustering_monte_carlo_standard_ci(
      data = sim_data,
      true_tau = as.numeric(scenario[["true_tau"]]),
      alpha = alpha
    )

    if (is.null(standard_result)) {
      next
    }

    valid_replications <- valid_replications + 1L
    if (identical(estimator, "cluster_robust_t_ci")) {
      metric_hits[[valid_replications]] <- isTRUE(
        standard_result$covers
      )
    } else {
      t_stat <- standard_result$att / standard_result$se
      p_value <- 2 * stats::pt(-abs(t_stat), df = as.integer(scenario[["G"]]) - 1L)
      metric_hits[[valid_replications]] <- isTRUE(p_value < alpha)
    }
  }

  metric_name <- switch(
    estimator,
    cluster_robust_t_ci = "coverage_rate",
    cluster_robust_t_test = if (isTRUE(all.equal(
      as.numeric(scenario[["true_tau"]]),
      0
    ))) {
      "rejection_rate"
    } else {
      "power"
    }
  )

  list(
    scenario_id = scenario[["id"]],
    n_simulations = n_target,
    n_attempts = attempts,
    metric_name = metric_name,
    metric_value = mean(metric_hits)
  )
}
