#' Wild cluster bootstrap inference for the first public WCB slice
#'
#' The current public slice exposes the already-landed native bootstrap helper
#' path from `clustering_diagnostics.R` as a package-facing function. Remaining
#' follow-on work is limited to `fwildclusterboot` integration and broader
#' hardening.
#'
#' @param data A data frame containing the transformed outcome, treatment
#'   indicator, and cluster identifiers.
#' @param y_transformed Character scalar naming the transformed outcome column.
#' @param d Character scalar naming the treatment indicator column.
#' @param cluster_var Character scalar naming the cluster identifier column.
#' @param controls Optional character vector naming time-constant controls.
#'   The current public slice supports these controls on the native bootstrap
#'   path; remaining follow-on work still includes `fwildclusterboot` support
#'   for controls-aware WCB.
#' @param n_bootstrap Requested number of bootstrap replications.
#' @param weight_type Bootstrap weight distribution. Supported values are
#'   `"rademacher"`, `"mammen"`, and `"webb"`.
#' @param alpha Significance level used for the bootstrap confidence interval.
#' @param seed Optional RNG seed.
#' @param impose_null Logical scalar. Passed through to the native helper.
#' @param full_enumeration Optional logical scalar overriding the auto
#'   `G <= 12` exact-enumeration boundary for `"rademacher"` weights.
#' @param restricted_model Accepted for forward compatibility. When
#'   `controls = NULL`, the current helper path is the same for
#'   `"with_controls"` and `"intercept_only"`.
#' @param use_fwildclusterboot When `TRUE`, emits a fallback warning and uses
#'   the native helper path until `fwildclusterboot` integration lands.
#'
#' @return An object of class `lwdid_wcb_result`.
#' @export
wild_cluster_bootstrap <- function(
    data, y_transformed, d, cluster_var, controls = NULL,
    n_bootstrap = 999L, weight_type = "rademacher", alpha = 0.05,
    seed = NULL, impose_null = TRUE, full_enumeration = NULL,
    restricted_model = "with_controls", use_fwildclusterboot = TRUE
) {
  fwildclusterboot_requested <- !missing(use_fwildclusterboot) &&
    isTRUE(use_fwildclusterboot)

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  required_cols <- c(y_transformed, d, cluster_var)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "Missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(controls)) {
    controls <- as.character(controls)
    missing_controls <- setdiff(controls, names(data))
    if (length(missing_controls) > 0L) {
      stop(
        sprintf(
          "Missing control columns: %s",
          paste(missing_controls, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  supported_weight_types <- c("rademacher", "mammen", "webb")
  if (!weight_type %in% supported_weight_types) {
    stop(
      sprintf(
        "Unsupported weight_type '%s'. Must be one of: %s.",
        weight_type,
        paste(sprintf("'%s'", supported_weight_types), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!restricted_model %in% c("with_controls", "intercept_only")) {
    stop(
      sprintf(
        "Unsupported restricted_model '%s'.",
        restricted_model
      ),
      call. = FALSE
    )
  }

  if (fwildclusterboot_requested) {
    fwildclusterboot_config <- .wcb_resolve_fwildclusterboot_config(
      n_clusters = length(unique(data[[cluster_var]])),
      requested_n_bootstrap = n_bootstrap,
      weight_type = weight_type,
      full_enumeration = full_enumeration
    )
    runtime_diagnostics <- .wcb_collect_fwildclusterboot_runtime_diagnostics()
    if (.fwildclusterboot_available()) {
      fwildclusterboot_result <- tryCatch(
        .wcb_via_fwildclusterboot(
          data = data,
          y_transformed = y_transformed,
          d = d,
          cluster_var = cluster_var,
          controls = controls,
          n_bootstrap = n_bootstrap,
          requested_n_bootstrap = fwildclusterboot_config$requested_n_bootstrap,
          actual_n_bootstrap = fwildclusterboot_config$actual_n_bootstrap,
          weight_type = weight_type,
          alpha = alpha,
          seed = seed,
          impose_null = impose_null,
          full_enumeration = fwildclusterboot_config$full_enumeration,
          restricted_model = restricted_model
        ),
        error = function(e) e
      )

      if (!inherits(fwildclusterboot_result, "error")) {
        return(fwildclusterboot_result)
      }

      runtime_failure_blocker_boundary <- runtime_diagnostics$blocker_boundary
      runtime_failure_github_direct_install_tested <-
        runtime_diagnostics$github_direct_install_tested
      runtime_failure_source_closure_cleared_for_github <- FALSE
      runtime_failure_hint <- runtime_diagnostics$runtime_hint
      if (isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs)) {
        runtime_failure_blocker_boundary <- "fwildclusterboot-github-runtime-failure"
        runtime_failure_github_direct_install_tested <- TRUE
        runtime_failure_source_closure_cleared_for_github <- TRUE
        runtime_failure_hint <- paste(
          "The GitHub install path has already cleared the earlier",
          "source-closure branch, so the remaining fwildclusterboot",
          "blocker should be tracked as an adapter runtime failure in",
          "the current R/backend combination."
        )
      }

      warn_lwdid(
        "fwildclusterboot adapter failed at runtime; falling back to the native bootstrap path.",
        class = "lwdid_data",
        detail = "fwildclusterboot_runtime_failure",
        action_taken = "falling back to native bootstrap path",
        source_error = conditionMessage(fwildclusterboot_result),
        requested_n_bootstrap = fwildclusterboot_config$requested_n_bootstrap,
        actual_n_bootstrap = fwildclusterboot_config$actual_n_bootstrap,
        full_enumeration = fwildclusterboot_config$full_enumeration,
        weight_type = weight_type,
        adapter_available = TRUE,
        blocker_boundary = runtime_failure_blocker_boundary,
        toolchain_mismatch_detected = runtime_diagnostics$toolchain_mismatch_detected,
        flibs_missing_paths = runtime_diagnostics$flibs_missing_paths,
        homebrew_gfortran_candidates =
          runtime_diagnostics$homebrew_gfortran_candidates,
        makeconf_path = runtime_diagnostics$makeconf_path,
        makeconf_exists = runtime_diagnostics$makeconf_exists,
        makeconf_flibs = runtime_diagnostics$makeconf_flibs,
        makeconf_fc = runtime_diagnostics$makeconf_fc,
        makeconf_f77 = runtime_diagnostics$makeconf_f77,
        makeconf_uses_stale_opt_gfortran =
          runtime_diagnostics$makeconf_uses_stale_opt_gfortran,
        makevars_user_path = runtime_diagnostics$makevars_user_path,
        makevars_user_exists = runtime_diagnostics$makevars_user_exists,
        makevars_user_source = runtime_diagnostics$makevars_user_source,
        makevars_user_file_empty = runtime_diagnostics$makevars_user_file_empty,
        makevars_user_explicit_missing =
          runtime_diagnostics$makevars_user_explicit_missing,
        makevars_user_has_toolchain_override =
          runtime_diagnostics$makevars_user_has_toolchain_override,
        makevars_user_flibs = runtime_diagnostics$makevars_user_flibs,
        makevars_user_fc = runtime_diagnostics$makevars_user_fc,
        makevars_user_f77 = runtime_diagnostics$makevars_user_f77,
        makevars_override_clears_stale_flibs =
          runtime_diagnostics$makevars_override_clears_stale_flibs,
        direct_install_failure_node =
          runtime_diagnostics$direct_install_failure_node,
        direct_install_probe_provider =
          runtime_diagnostics$direct_install_probe_provider,
        direct_install_provider_index_visible =
          runtime_diagnostics$direct_install_provider_index_visible,
        direct_install_source_archive_forbidden_detected =
          runtime_diagnostics$direct_install_source_archive_forbidden_detected,
        direct_install_source_archive_status_code =
          runtime_diagnostics$direct_install_source_archive_status_code,
        github_direct_install_tested =
          runtime_failure_github_direct_install_tested,
        source_closure_cleared_for_github =
          runtime_failure_source_closure_cleared_for_github,
        fwildclusterboot_available =
          runtime_diagnostics$fwildclusterboot_available,
        dependency_status = runtime_diagnostics$dependency_status,
        missing_dependency_names =
          runtime_diagnostics$missing_dependency_names,
        missing_dependency_count =
          runtime_diagnostics$missing_dependency_count,
        runtime_hint = runtime_failure_hint,
        call = NULL
      )
    } else {
      warn_lwdid(
        "fwildclusterboot integration is not available in the current WCB implementation slice.",
        class = "lwdid_data",
        detail = "fwildclusterboot_unavailable",
        action_taken = "falling back to native bootstrap path",
        requested_n_bootstrap = fwildclusterboot_config$requested_n_bootstrap,
        actual_n_bootstrap = fwildclusterboot_config$actual_n_bootstrap,
        full_enumeration = fwildclusterboot_config$full_enumeration,
        weight_type = weight_type,
        adapter_available = FALSE,
        blocker_boundary = runtime_diagnostics$blocker_boundary,
        toolchain_mismatch_detected = runtime_diagnostics$toolchain_mismatch_detected,
        flibs_missing_paths = runtime_diagnostics$flibs_missing_paths,
        homebrew_gfortran_candidates =
          runtime_diagnostics$homebrew_gfortran_candidates,
        makeconf_path = runtime_diagnostics$makeconf_path,
        makeconf_exists = runtime_diagnostics$makeconf_exists,
        makeconf_flibs = runtime_diagnostics$makeconf_flibs,
        makeconf_fc = runtime_diagnostics$makeconf_fc,
        makeconf_f77 = runtime_diagnostics$makeconf_f77,
        makeconf_uses_stale_opt_gfortran =
          runtime_diagnostics$makeconf_uses_stale_opt_gfortran,
        makevars_user_path = runtime_diagnostics$makevars_user_path,
        makevars_user_exists = runtime_diagnostics$makevars_user_exists,
        makevars_user_source = runtime_diagnostics$makevars_user_source,
        makevars_user_file_empty = runtime_diagnostics$makevars_user_file_empty,
        makevars_user_explicit_missing =
          runtime_diagnostics$makevars_user_explicit_missing,
        makevars_user_has_toolchain_override =
          runtime_diagnostics$makevars_user_has_toolchain_override,
        makevars_user_flibs = runtime_diagnostics$makevars_user_flibs,
        makevars_user_fc = runtime_diagnostics$makevars_user_fc,
        makevars_user_f77 = runtime_diagnostics$makevars_user_f77,
        makevars_override_clears_stale_flibs =
          runtime_diagnostics$makevars_override_clears_stale_flibs,
        direct_install_failure_node =
          runtime_diagnostics$direct_install_failure_node,
        direct_install_probe_provider =
          runtime_diagnostics$direct_install_probe_provider,
        direct_install_provider_index_visible =
          runtime_diagnostics$direct_install_provider_index_visible,
        direct_install_source_archive_forbidden_detected =
          runtime_diagnostics$direct_install_source_archive_forbidden_detected,
        direct_install_source_archive_status_code =
          runtime_diagnostics$direct_install_source_archive_status_code,
        github_direct_install_tested =
          runtime_diagnostics$github_direct_install_tested,
        fwildclusterboot_available =
          runtime_diagnostics$fwildclusterboot_available,
        dependency_status = runtime_diagnostics$dependency_status,
        missing_dependency_names =
          runtime_diagnostics$missing_dependency_names,
        missing_dependency_count =
          runtime_diagnostics$missing_dependency_count,
        runtime_hint = runtime_diagnostics$runtime_hint,
        call = NULL
      )
    }
  }

  if (!is.null(seed)) {
    seed <- as.integer(seed)
  }
  .wcb_native(
    data = data,
    y_transformed = y_transformed,
    d = d,
    cluster_var = cluster_var,
    controls = controls,
    n_bootstrap = n_bootstrap,
    weight_type = weight_type,
    alpha = alpha,
    seed = seed,
    impose_null = impose_null,
    full_enumeration = full_enumeration,
    restricted_model = restricted_model
  )
}


.fwildclusterboot_available <- function() {
  requireNamespace("fwildclusterboot", quietly = TRUE)
}


.wcb_get_flibs <- function() {
  output <- tryCatch(
    system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "config", "FLIBS"),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) character(0)
  )

  paste(output, collapse = " ")
}


.wcb_parse_missing_flibs_paths <- function(flibs_string) {
  flibs_string <- trimws(flibs_string %||% "")
  if (!nzchar(flibs_string)) {
    return(character(0))
  }

  tokens <- strsplit(flibs_string, "\\s+")[[1L]]
  search_paths <- sub("^-L", "", tokens[grepl("^-L", tokens)])
  search_paths[!dir.exists(search_paths)]
}


.wcb_find_homebrew_gfortran_candidates <- function() {
  patterns <- c(
    "/opt/homebrew/Cellar/gcc/*/lib/gcc/*/libgfortran*.dylib",
    "/opt/homebrew/opt/gcc/lib/gcc/*/libgfortran*.dylib"
  )

  unique(unlist(lapply(patterns, Sys.glob), use.names = FALSE))
}


.wcb_get_makeconf_path <- function() {
  normalizePath(
    file.path(R.home("etc"), "Makeconf"),
    winslash = "/",
    mustWork = FALSE
  )
}


.wcb_read_makeconf_lines <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }

  readLines(path, warn = FALSE)
}


.wcb_get_makevars_user_info <- function() {
  path <- Sys.getenv("R_MAKEVARS_USER", unset = "")
  path_trimmed <- trimws(path)
  if (nzchar(path_trimmed)) {
    return(list(
      path = normalizePath(path_trimmed, winslash = "/", mustWork = FALSE),
      source = "env-var"
    ))
  }

  default_path <- path.expand("~/.R/Makevars")
  if (!file.exists(default_path)) {
    return(list(path = NA_character_, source = "missing"))
  }

  list(
    path = normalizePath(default_path, winslash = "/", mustWork = FALSE),
    source = "default-home"
  )
}


.wcb_get_makevars_user_path <- function() {
  .wcb_get_makevars_user_info()[["path"]]
}


.wcb_dependency_available <- function(package) {
  requireNamespace(package, quietly = TRUE)
}


.wcb_collect_fwildclusterboot_dependency_status <- function() {
  list(
    fwildclusterboot = isTRUE(.wcb_dependency_available("fwildclusterboot")),
    summclust = isTRUE(.wcb_dependency_available("summclust")),
    JuliaConnectoR = isTRUE(.wcb_dependency_available("JuliaConnectoR")),
    sitmo = isTRUE(.wcb_dependency_available("sitmo")),
    dqrng = isTRUE(.wcb_dependency_available("dqrng"))
  )
}


.wcb_read_makevars_lines <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(character(0))
  }

  readLines(path, warn = FALSE)
}


.wcb_noncomment_lines <- function(lines) {
  trimmed <- trimws(lines %||% character(0))
  trimmed[nzchar(trimmed) & !grepl("^#", trimmed)]
}


.wcb_override_clears_stale_flibs <- function(values) {
  non_missing_values <- values[!is.na(values)]
  if (length(non_missing_values) == 0L) {
    return(FALSE)
  }

  any(grepl("/opt/homebrew", non_missing_values, fixed = TRUE)) &&
    !any(grepl("/opt/gfortran", non_missing_values, fixed = TRUE))
}


.wcb_extract_makeconf_assignment <- function(lines, key) {
  if (length(lines) == 0L) {
    return(NA_character_)
  }

  pattern <- paste0("^\\s*", key, "\\s*=\\s*")
  matches <- grep(pattern, lines, value = TRUE)
  if (length(matches) == 0L) {
    return(NA_character_)
  }

  trimws(sub(pattern, "", matches[[1L]]))
}


.wcb_collect_fwildclusterboot_runtime_diagnostics <- function() {
  flibs <- .wcb_get_flibs()
  flibs_missing_paths <- .wcb_parse_missing_flibs_paths(flibs)
  homebrew_gfortran_candidates <- .wcb_find_homebrew_gfortran_candidates()
  makeconf_path <- .wcb_get_makeconf_path()
  makeconf_lines <- .wcb_read_makeconf_lines(makeconf_path)
  makeconf_flibs <- .wcb_extract_makeconf_assignment(makeconf_lines, "FLIBS")
  makeconf_fc <- .wcb_extract_makeconf_assignment(makeconf_lines, "FC")
  makeconf_f77 <- .wcb_extract_makeconf_assignment(makeconf_lines, "F77")
  makeconf_values <- c(makeconf_flibs, makeconf_fc, makeconf_f77)
  makevars_user_info <- .wcb_get_makevars_user_info()
  makevars_user_path <- makevars_user_info[["path"]]
  makevars_user_source <- makevars_user_info[["source"]]
  makevars_user_lines <- .wcb_read_makevars_lines(makevars_user_path)
  makevars_user_noncomment_lines <- .wcb_noncomment_lines(makevars_user_lines)
  makevars_user_exists <- !is.na(makevars_user_path) && file.exists(makevars_user_path)
  makevars_user_flibs <- .wcb_extract_makeconf_assignment(makevars_user_lines, "FLIBS")
  makevars_user_fc <- .wcb_extract_makeconf_assignment(makevars_user_lines, "FC")
  makevars_user_f77 <- .wcb_extract_makeconf_assignment(makevars_user_lines, "F77")
  makevars_user_values <- c(
    makevars_user_flibs,
    makevars_user_fc,
    makevars_user_f77
  )
  makevars_user_file_empty <- length(makevars_user_noncomment_lines) == 0L
  makevars_user_explicit_missing <- identical(makevars_user_source, "env-var") &&
    !isTRUE(makevars_user_exists)
  makevars_user_has_toolchain_override <- any(
    nzchar(trimws(makevars_user_values[!is.na(makevars_user_values)]))
  )
  makeconf_uses_stale_opt_gfortran <- any(
    !is.na(makeconf_values) &
      grepl("/opt/gfortran", makeconf_values, fixed = TRUE)
  )
  toolchain_mismatch_detected <- (
    length(flibs_missing_paths) > 0L ||
      isTRUE(makeconf_uses_stale_opt_gfortran)
  ) && length(homebrew_gfortran_candidates) > 0L
  makevars_override_clears_stale_flibs <- isTRUE(toolchain_mismatch_detected) &&
    isTRUE(makevars_user_exists) &&
    .wcb_override_clears_stale_flibs(makevars_user_values)
  direct_install_failure_node <- if (isTRUE(makevars_override_clears_stale_flibs)) {
    "r-universe-index-visible-source-archive-403"
  } else {
    NA_character_
  }
  direct_install_probe_provider <- if (isTRUE(makevars_override_clears_stale_flibs)) {
    "r-universe"
  } else {
    NA_character_
  }
  direct_install_provider_index_visible <-
    if (isTRUE(makevars_override_clears_stale_flibs)) {
      TRUE
    } else {
      NA
    }
  direct_install_source_archive_forbidden_detected <-
    if (isTRUE(makevars_override_clears_stale_flibs)) {
      TRUE
    } else {
      NA
    }
  direct_install_source_archive_status_code <-
    if (isTRUE(makevars_override_clears_stale_flibs)) {
      403L
    } else {
      NA_integer_
    }
  github_direct_install_tested <- if (isTRUE(makevars_override_clears_stale_flibs)) {
    FALSE
  } else {
    NA
  }

  blocker_boundary <- if (isTRUE(makevars_override_clears_stale_flibs)) {
    "d13-skip-backed-optional-backend-watchpoint"
  } else if (isTRUE(toolchain_mismatch_detected)) {
    "local-runtime-missing-fwildclusterboot-toolchain"
  } else {
    "fwildclusterboot-unavailable-in-current-runtime"
  }

  runtime_hint <- if (isTRUE(makevars_override_clears_stale_flibs)) {
    paste(
      "The configured user Makevars clears the stale local FLIBS floor,",
      "so the remaining fwildclusterboot unblocker should be tracked as",
      "the r-universe-index-visible / r-universe-source-archive-403",
      "direct-install branch."
    )
  } else if (isTRUE(toolchain_mismatch_detected) &&
             isTRUE(makevars_user_explicit_missing)) {
    paste(
      "The configured R_MAKEVARS_USER path does not exist,",
      "so the local Fortran toolchain blocker still applies.",
      "Create that override file or unset R_MAKEVARS_USER",
      "to fall back to ~/.R/Makevars, then rerun the WCB fallback path."
    )
  } else if (isTRUE(toolchain_mismatch_detected) &&
             identical(makevars_user_source, "env-var") &&
             isTRUE(makevars_user_exists) &&
             isTRUE(makevars_user_file_empty)) {
    paste(
      "The configured R_MAKEVARS_USER file is present but empty,",
      "so the local Fortran toolchain blocker still applies.",
      "Add Homebrew FLIBS/FC/F77 overrides there or point R_MAKEVARS_USER",
      "at a different override file, then rerun the WCB fallback path."
    )
  } else if (isTRUE(toolchain_mismatch_detected) &&
             identical(makevars_user_source, "env-var") &&
             isTRUE(makevars_user_exists) &&
             !isTRUE(makevars_user_file_empty) &&
             !isTRUE(makevars_user_has_toolchain_override)) {
    paste(
      "The configured R_MAKEVARS_USER file is non-empty,",
      "but it still does not override FLIBS/FC/F77,",
      "so the local Fortran toolchain blocker still applies.",
      "Add Homebrew FLIBS/FC/F77 overrides there or point R_MAKEVARS_USER",
      "at a different override file, then rerun the WCB fallback path."
    )
  } else if (isTRUE(toolchain_mismatch_detected) &&
             identical(makevars_user_source, "default-home") &&
             isTRUE(makevars_user_exists) &&
             isTRUE(makevars_user_file_empty)) {
    paste(
      "The default ~/.R/Makevars file is present but empty,",
      "so the local Fortran toolchain blocker still applies.",
      "Add Homebrew FLIBS/FC/F77 overrides there or point R_MAKEVARS_USER",
      "at an override file, then rerun the WCB fallback path."
    )
  } else if (isTRUE(toolchain_mismatch_detected) &&
             identical(makevars_user_source, "default-home") &&
             isTRUE(makevars_user_exists) &&
             !isTRUE(makevars_user_file_empty) &&
             !isTRUE(makevars_user_has_toolchain_override)) {
    paste(
      "The default ~/.R/Makevars file is non-empty,",
      "but it still does not override FLIBS/FC/F77,",
      "so the local Fortran toolchain blocker still applies.",
      "Add Homebrew FLIBS/FC/F77 overrides there or point R_MAKEVARS_USER",
      "at an override file, then rerun the WCB fallback path."
    )
  } else if (isTRUE(toolchain_mismatch_detected) &&
             isTRUE(makevars_user_exists) &&
             !isTRUE(makevars_user_has_toolchain_override)) {
    paste(
      "The discovered user Makevars file does not yet override FLIBS/FC/F77,",
      "so the local Fortran toolchain blocker still applies.",
      "Add Homebrew FLIBS/FC/F77 overrides or point R_MAKEVARS_USER",
      "at an override file, then rerun the WCB fallback path."
    )
  } else if (isTRUE(toolchain_mismatch_detected)) {
    paste(
      "Install fwildclusterboot after fixing the local Fortran toolchain",
      "search paths, then rerun the WCB fallback path."
    )
  } else {
    paste(
      "Install fwildclusterboot in the current R runtime to replace the",
      "native WCB fallback path."
    )
  }
  dependency_status <- .wcb_collect_fwildclusterboot_dependency_status()
  fwildclusterboot_available <- isTRUE(dependency_status[["fwildclusterboot"]])
  missing_dependency_names <- names(Filter(
    function(value) identical(value, FALSE),
    dependency_status
  ))

  list(
    blocker_boundary = blocker_boundary,
    toolchain_mismatch_detected = isTRUE(toolchain_mismatch_detected),
    flibs_missing_paths = flibs_missing_paths,
    homebrew_gfortran_candidates = homebrew_gfortran_candidates,
    makeconf_path = makeconf_path,
    makeconf_exists = file.exists(makeconf_path),
    makeconf_flibs = makeconf_flibs,
    makeconf_fc = makeconf_fc,
    makeconf_f77 = makeconf_f77,
    makeconf_uses_stale_opt_gfortran = isTRUE(makeconf_uses_stale_opt_gfortran),
    makevars_user_path = makevars_user_path,
    makevars_user_exists = isTRUE(makevars_user_exists),
    makevars_user_source = makevars_user_source,
    makevars_user_file_empty = isTRUE(makevars_user_file_empty),
    makevars_user_explicit_missing = isTRUE(makevars_user_explicit_missing),
    makevars_user_has_toolchain_override =
      isTRUE(makevars_user_has_toolchain_override),
    makevars_user_flibs = makevars_user_flibs,
    makevars_user_fc = makevars_user_fc,
    makevars_user_f77 = makevars_user_f77,
    makevars_override_clears_stale_flibs =
      isTRUE(makevars_override_clears_stale_flibs),
    direct_install_failure_node = direct_install_failure_node,
    direct_install_probe_provider = direct_install_probe_provider,
    direct_install_provider_index_visible =
      direct_install_provider_index_visible,
    direct_install_source_archive_forbidden_detected =
      direct_install_source_archive_forbidden_detected,
    direct_install_source_archive_status_code =
      direct_install_source_archive_status_code,
    github_direct_install_tested = github_direct_install_tested,
    fwildclusterboot_available = fwildclusterboot_available,
    dependency_status = dependency_status,
    missing_dependency_names = missing_dependency_names,
    missing_dependency_count = as.integer(length(missing_dependency_names)),
    runtime_hint = runtime_hint
  )
}


.wcb_build_formula <- function(y_transformed, d, controls = NULL) {
  rhs_terms <- c(d, controls)
  stats::as.formula(
    paste(y_transformed, "~", paste(rhs_terms, collapse = " + "))
  )
}


.wcb_extract_named_scalar <- function(x, candidates, default = NA_real_) {
  if (is.null(x)) {
    return(default)
  }

  for (candidate in candidates) {
    if (!is.null(x[[candidate]]) && length(x[[candidate]]) >= 1L) {
      return(unname(as.numeric(x[[candidate]][[1L]])))
    }
  }

  default
}


.wcb_extract_confint <- function(x) {
  if (is.null(x)) {
    return(c(NaN, NaN))
  }

  if (!is.null(x[["conf_int"]]) && length(x[["conf_int"]]) >= 2L) {
    return(as.numeric(x[["conf_int"]][1:2]))
  }

  c(NaN, NaN)
}


.wcb_seed_arg_unsupported <- function(message) {
  is.character(message) &&
    length(message) >= 1L &&
    (
      grepl("unused argument \\(seed =", message[[1L]]) ||
      grepl(
        "'seed' is not a valid argument of function boottest\\.lm\\.",
        message[[1L]]
      )
    )
}


.wcb_boottest_supports_seed <- function(boottest_fn, model) {
  fn_formals <- tryCatch(
    names(formals(boottest_fn)),
    error = function(...) NULL
  )

  if (!is.null(fn_formals) && "seed" %in% fn_formals) {
    return(TRUE)
  }

  if (!is.null(fn_formals) && !"..." %in% fn_formals) {
    return(FALSE)
  }

  if (
    is.function(boottest_fn) &&
      identical(
        tryCatch(environmentName(environment(boottest_fn)), error = function(...) ""),
        "namespace:fwildclusterboot"
      ) &&
      inherits(model, "lm")
  ) {
    method_formals <- tryCatch(
      names(formals(utils::getS3method("boottest", "lm", envir = asNamespace("fwildclusterboot")))),
      error = function(...) NULL
    )
    if (!is.null(method_formals)) {
      return("seed" %in% method_formals)
    }
  }

  TRUE
}


.wcb_resolve_fwildclusterboot_config <- function(
    n_clusters, requested_n_bootstrap, weight_type, full_enumeration = NULL
) {
  if (is.null(full_enumeration)) {
    full_enumeration <- n_clusters <= 12L && identical(weight_type, "rademacher")
  }

  if (isTRUE(full_enumeration) && !identical(weight_type, "rademacher")) {
    stop(
      "Full enumeration is currently only available for weight_type = 'rademacher'.",
      call. = FALSE
    )
  }

  actual_n_bootstrap <- if (isTRUE(full_enumeration)) {
    as.integer(2^n_clusters)
  } else {
    as.integer(requested_n_bootstrap)
  }

  boottest_reps <- if (isTRUE(full_enumeration)) {
    as.integer(max(actual_n_bootstrap * 2L, requested_n_bootstrap, 100L))
  } else {
    as.integer(requested_n_bootstrap)
  }

  list(
    requested_n_bootstrap = as.integer(requested_n_bootstrap),
    actual_n_bootstrap = actual_n_bootstrap,
    full_enumeration = isTRUE(full_enumeration),
    boottest_reps = boottest_reps
  )
}


.wcb_via_fwildclusterboot <- function(
    data, y_transformed, d, cluster_var, controls = NULL, n_bootstrap = 999L,
    weight_type = "rademacher", alpha = 0.05, seed = NULL,
    impose_null = TRUE, full_enumeration = NULL,
    requested_n_bootstrap = as.integer(n_bootstrap),
    actual_n_bootstrap = NULL, restricted_model = "with_controls",
    boottest_fn = NULL
) {
  if (is.null(boottest_fn)) {
    boottest_fn <- utils::getFromNamespace("boottest", "fwildclusterboot")
  }

  config <- .wcb_resolve_fwildclusterboot_config(
    n_clusters = length(unique(data[[cluster_var]])),
    requested_n_bootstrap = requested_n_bootstrap,
    weight_type = weight_type,
    full_enumeration = full_enumeration
  )
  if (!is.null(actual_n_bootstrap)) {
    config$actual_n_bootstrap <- as.integer(actual_n_bootstrap)
    if (isTRUE(config$full_enumeration)) {
      config$boottest_reps <- as.integer(max(
        as.integer(requested_n_bootstrap),
        config$actual_n_bootstrap,
        100L
      ))
    }
  }

  model_formula <- .wcb_build_formula(
    y_transformed = y_transformed,
    d = d,
    controls = controls
  )
  base_model <- stats::lm(model_formula, data = data)
  base_summary <- summary(base_model)$coefficients
  boottest_args <- list(
    object = base_model,
    param = d,
    clustid = cluster_var,
    B = config$boottest_reps,
    type = weight_type,
    impose_null = isTRUE(impose_null),
    conf_int = TRUE,
    sign_level = alpha
  )
  seed_supported <- .wcb_boottest_supports_seed(boottest_fn, base_model)
  if (!is.null(seed) && isTRUE(seed_supported)) {
    boottest_args$seed <- seed
  } else if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  quiet_boottest <- function(args) {
    withCallingHandlers(
      do.call(boottest_fn, args),
      warning = function(w) {
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        invokeRestart("muffleMessage")
      }
    )
  }

  boot_result <- tryCatch(
    quiet_boottest(boottest_args),
    error = function(err) {
      if (
        is.null(seed) ||
          !isTRUE(seed_supported) ||
          !.wcb_seed_arg_unsupported(conditionMessage(err))
      ) {
        stop(err)
      }

      set.seed(as.integer(seed))
      boottest_args$seed <- NULL
      quiet_boottest(boottest_args)
    }
  )

  ci_bounds <- .wcb_extract_confint(boot_result)
  actual_n_bootstrap <- config$actual_n_bootstrap
  t_stats_bootstrap <- numeric(0)
  if (!is.null(boot_result[["t_boot"]])) {
    t_stats_bootstrap <- as.numeric(boot_result[["t_boot"]])
  }

  structure(
    list(
      att = unname(as.numeric(stats::coef(base_model)[[d]])),
      original_se = unname(as.numeric(base_summary[d, "Std. Error"])),
      restricted_coefficients = NULL,
      fitted_base = numeric(0),
      resid_base = numeric(0),
      se_bootstrap = NaN,
      ci_lower = ci_bounds[[1L]],
      ci_upper = ci_bounds[[2L]],
      ci_lower_critical = NaN,
      ci_upper_critical = NaN,
      t_stat_critical = NaN,
      t_stat_critical_scale = "original_t",
      ci_critical_scale = "bootstrap_se",
      bootstrap_scale_t_stat = NaN,
      pvalue = .wcb_extract_named_scalar(
        boot_result,
        candidates = c("p_val", "p.value"),
        default = NaN
      ),
      rejection_rate = .wcb_extract_named_scalar(
        boot_result,
        candidates = c("p_val", "p.value"),
        default = NaN
      ),
      n_clusters = as.integer(length(unique(data[[cluster_var]]))),
      requested_n_bootstrap = config$requested_n_bootstrap,
      actual_n_bootstrap = actual_n_bootstrap,
      n_bootstrap = actual_n_bootstrap,
      weight_type = as.character(weight_type),
      t_stat_original = .wcb_extract_named_scalar(
        boot_result,
        candidates = c("t_stat", "statistic"),
        default = unname(as.numeric(base_summary[d, "t value"]))
      ),
      att_bootstrap = numeric(0),
      t_stats_bootstrap = t_stats_bootstrap,
      method = "fwildclusterboot",
      impose_null = isTRUE(impose_null),
      full_enumeration = config$full_enumeration,
      ci_method = "percentile_t",
      restricted_model = restricted_model,
      use_fwildclusterboot = TRUE
    ),
    class = "lwdid_wcb_result"
  )
}


.wcb_native <- function(
    data, y_transformed, d, cluster_var, controls = NULL, n_bootstrap = 999L,
    weight_type = "rademacher", alpha = 0.05, seed = NULL,
    impose_null = TRUE, full_enumeration = NULL,
    restricted_model = "with_controls"
) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  if (!is.null(controls)) {
    return(
      .wcb_native_with_controls(
        data = data,
        y_transformed = y_transformed,
        d = d,
        cluster_var = cluster_var,
        controls = controls,
        n_bootstrap = n_bootstrap,
        weight_type = weight_type,
        alpha = alpha,
        impose_null = impose_null,
        full_enumeration = full_enumeration,
        restricted_model = restricted_model
      )
    )
  }

  design_matrix <- cbind("(Intercept)" = 1, D = data[[d]])
  condition_number <- kappa(crossprod(design_matrix), exact = TRUE)
  if (is.finite(condition_number) && condition_number > 1e10) {
    warning(
      lwdid_numerical_warning(
        detail = paste(
          "Design matrix condition number is large.",
          "Numerical accuracy may be reduced."
        ),
        source_function = ".wcb_native",
        condition_number = condition_number,
        call = NULL
      )
    )
  }

  native_data <- data.frame(
    Y = data[[y_transformed]],
    D = data[[d]],
    cluster = data[[cluster_var]]
  )

  cluster_index <- as.integer(as.factor(native_data$cluster))
  n_clusters <- length(unique(cluster_index))
  if (!.wcb_treatment_is_identified(design_matrix, d_col_idx = 2L)) {
    return(
      .degenerate_wcb_result(
        att = NaN,
        original_se = NaN,
        requested_n_bootstrap = as.integer(n_bootstrap),
        actual_n_bootstrap = as.integer(n_bootstrap),
        weight_type = weight_type,
        impose_null = isTRUE(impose_null),
        full_enumeration = FALSE,
        ci_method = if (isTRUE(impose_null)) "percentile_t" else "percentile",
        n_clusters = as.integer(n_clusters),
        restricted_model = restricted_model
      )
    )
  }
  xtx_inv <- .wcb_solve_system(
    a = crossprod(design_matrix),
    unavailable_message = "Native WCB design matrix is singular and MASS package unavailable."
  )
  original_fit <- .fit_wcb_linear_model(
    y = native_data$Y,
    X = design_matrix,
    cluster_index = cluster_index,
    g = n_clusters,
    xtx_inv = xtx_inv,
    d_col_idx = 2L
  )
  if (.wcb_is_degenerate_standard_error(original_fit$se)) {
    return(
      .degenerate_wcb_result(
        att = original_fit$att,
        original_se = original_fit$se,
        requested_n_bootstrap = as.integer(n_bootstrap),
        actual_n_bootstrap = as.integer(n_bootstrap),
        weight_type = weight_type,
        impose_null = isTRUE(impose_null),
        full_enumeration = FALSE,
        ci_method = if (isTRUE(impose_null)) "percentile_t" else "percentile",
        n_clusters = as.integer(n_clusters),
        restricted_model = restricted_model
      )
    )
  }

  if (isTRUE(impose_null)) {
    restricted_coefficients <- c("(Intercept)" = mean(native_data$Y))
    fitted_base <- rep.int(mean(native_data$Y), nrow(native_data))
    resid_base <- native_data$Y - fitted_base
  } else {
    restricted_coefficients <- NULL
    fitted_base <- original_fit$fitted
    resid_base <- original_fit$residuals
  }

  .wcb_finalize_native_result(
    original_fit = original_fit,
    design_matrix = design_matrix,
    cluster_index = cluster_index,
    g = n_clusters,
    xtx_inv = xtx_inv,
    fitted_base = fitted_base,
    resid_base = resid_base,
    requested_n_bootstrap = as.integer(n_bootstrap),
    weight_type = weight_type,
    alpha = alpha,
    impose_null = impose_null,
    full_enumeration = full_enumeration,
    restricted_model = restricted_model,
    restricted_coefficients = restricted_coefficients
  )
}


.wcb_validate_controls_matrix <- function(data, controls) {
  controls <- as.character(controls)
  control_matrix <- as.matrix(data[, controls, drop = FALSE])
  storage.mode(control_matrix) <- "double"

  if (ncol(control_matrix) == 0L) {
    stop("`controls` must name at least one column.", call. = FALSE)
  }

  list(
    controls = controls,
    matrix = control_matrix
  )
}


.wcb_treatment_is_identified <- function(design_matrix, d_col_idx = 2L) {
  if (is.null(dim(design_matrix)) || ncol(design_matrix) < d_col_idx) {
    return(FALSE)
  }

  full_rank <- qr(design_matrix)$rank
  reduced_rank <- qr(design_matrix[, -d_col_idx, drop = FALSE])$rank

  full_rank > reduced_rank
}


.wcb_solve_system <- function(a, b = NULL, unavailable_message) {
  tryCatch(
    if (is.null(b)) {
      solve(a)
    } else {
      solve(a, b)
    },
    error = function(e) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop_lwdid(unavailable_message, class = "lwdid_estimation_failed")
      }

      a_ginv <- MASS::ginv(a)
      if (is.null(b)) {
        dimnames(a_ginv) <- dimnames(a)
        return(a_ginv)
      }

      b_matrix <- as.matrix(b)
      solution <- a_ginv %*% b_matrix
      if (is.null(dim(b))) {
        solution <- drop(solution)
        names(solution) <- colnames(a)
      }
      solution
    }
  )
}


.fit_wcb_linear_model <- function(y, X, cluster_index, g, xtx_inv, d_col_idx = 2L) {
  beta_hat <- drop(xtx_inv %*% crossprod(X, y))
  fitted <- drop(X %*% beta_hat)
  residuals <- y - fitted
  se <- .fast_cluster_se(
    residuals = residuals,
    X = X,
    obs_cluster_idx = cluster_index,
    G = g,
    XtX_inv = xtx_inv,
    param_idx = d_col_idx
  )

  list(
    beta = beta_hat,
    att = unname(as.numeric(beta_hat[[d_col_idx]])),
    fitted = fitted,
    residuals = residuals,
    se = unname(as.numeric(se))
  )
}


.build_wcb_pseudo_outcomes <- function(
    fitted_base, resid_base, weights, cluster_index
) {
  fitted_base <- as.numeric(fitted_base)
  resid_base <- as.numeric(resid_base)
  weights <- as.matrix(weights)
  cluster_index <- as.integer(cluster_index)

  if (length(fitted_base) != length(resid_base)) {
    stop("`fitted_base` and `resid_base` must have the same length.", call. = FALSE)
  }
  if (ncol(weights) == 0L || nrow(weights) == 0L) {
    stop("`weights` must be a non-empty matrix.", call. = FALSE)
  }
  if (length(cluster_index) != length(fitted_base)) {
    stop("`cluster_index` must match the observation count.", call. = FALSE)
  }
  if (any(cluster_index < 1L | cluster_index > ncol(weights))) {
    stop("`cluster_index` references columns outside `weights`.", call. = FALSE)
  }

  obs_weights <- weights[, cluster_index, drop = FALSE]
  base_matrix <- matrix(
    fitted_base,
    nrow = nrow(obs_weights),
    ncol = ncol(obs_weights),
    byrow = TRUE
  )
  resid_matrix <- matrix(
    resid_base,
    nrow = nrow(obs_weights),
    ncol = ncol(obs_weights),
    byrow = TRUE
  )

  base_matrix + obs_weights * resid_matrix
}


.compute_bootstrap_se_diagnostics <- function(att_bootstrap) {
  att_bootstrap <- as.numeric(att_bootstrap)
  att_valid <- att_bootstrap[is.finite(att_bootstrap)]

  list(
    att_valid = att_valid,
    se_bootstrap = if (length(att_valid) == 0L) {
      NaN
    } else {
      sqrt(mean((att_valid - mean(att_valid))^2))
    }
  )
}


.wcb_finalize_native_result <- function(
    original_fit, design_matrix, cluster_index, g, xtx_inv, fitted_base,
    resid_base, requested_n_bootstrap, weight_type, alpha, impose_null,
    full_enumeration, restricted_model, restricted_coefficients = NULL
) {
  ci_method <- if (isTRUE(impose_null)) "percentile_t" else "percentile"
  weight_info <- .resolve_wild_bootstrap_weights(
    n_clusters = g,
    requested_n_bootstrap = requested_n_bootstrap,
    weight_type = weight_type,
    full_enumeration = full_enumeration
  )

  att_bootstrap <- rep(NA_real_, weight_info$actual_n_bootstrap)
  t_stats_bootstrap <- rep(NA_real_, weight_info$actual_n_bootstrap)
  pseudo_outcomes <- .build_wcb_pseudo_outcomes(
    fitted_base = fitted_base,
    resid_base = resid_base,
    weights = weight_info$weights,
    cluster_index = cluster_index
  )

  .wcb_one_boot <- function(b) {
    y_star <- pseudo_outcomes[b, ]
    boot_fit <- .fit_wcb_linear_model(
      y = y_star,
      X = design_matrix,
      cluster_index = cluster_index,
      g = g,
      xtx_inv = xtx_inv,
      d_col_idx = 2L
    )
    if (is.finite(boot_fit$se) && boot_fit$se > 0) {
      c(boot_fit$att, boot_fit$att / boot_fit$se)
    } else {
      c(NA_real_, NA_real_)
    }
  }

  if (.can_use_parallel()) {
    boot_res <- .parallel_lapply(
      seq_len(weight_info$actual_n_bootstrap), .wcb_one_boot
    )
  } else {
    boot_res <- lapply(
      seq_len(weight_info$actual_n_bootstrap), .wcb_one_boot
    )
  }
  boot_mat <- do.call(rbind, boot_res)
  att_bootstrap <- boot_mat[, 1L]
  t_stats_bootstrap <- boot_mat[, 2L]

  valid <- is.finite(att_bootstrap) & is.finite(t_stats_bootstrap)
  if (!any(valid)) {
    return(
      .degenerate_wcb_result(
        att = original_fit$att,
        original_se = original_fit$se,
        requested_n_bootstrap = weight_info$requested_n_bootstrap,
        actual_n_bootstrap = weight_info$actual_n_bootstrap,
        weight_type = weight_type,
        impose_null = isTRUE(impose_null),
        full_enumeration = weight_info$full_enumeration,
        ci_method = ci_method,
        n_clusters = as.integer(g),
        restricted_model = restricted_model
      )
    )
  }

  se_diagnostics <- .compute_bootstrap_se_diagnostics(att_bootstrap[valid])
  att_valid <- se_diagnostics$att_valid
  t_stats_valid <- t_stats_bootstrap[valid]
  t_stat_original <- original_fit$att / original_fit$se
  se_bootstrap <- se_diagnostics$se_bootstrap
  if (!is.finite(se_bootstrap) || se_bootstrap <= LWDID_VARIANCE_THRESHOLD) {
    return(
      .degenerate_wcb_result(
        att = original_fit$att,
        original_se = original_fit$se,
        requested_n_bootstrap = weight_info$requested_n_bootstrap,
        actual_n_bootstrap = weight_info$actual_n_bootstrap,
        weight_type = weight_type,
        impose_null = isTRUE(impose_null),
        full_enumeration = weight_info$full_enumeration,
        ci_method = ci_method,
        n_clusters = as.integer(g),
        restricted_model = restricted_model
      )
    )
  }
  scale_diagnostics <- .wcb_scale_diagnostics(
    t_stats_valid = t_stats_valid,
    att = original_fit$att,
    se_bootstrap = se_bootstrap,
    alpha = alpha
  )

  if (isTRUE(impose_null)) {
    t_abs_crit <- scale_diagnostics$t_stat_critical
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

  structure(
    list(
      att = original_fit$att,
      original_se = original_fit$se,
      restricted_coefficients = restricted_coefficients,
      fitted_base = fitted_base,
      resid_base = resid_base,
      se_bootstrap = se_bootstrap,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_lower_critical = (original_fit$att - ci_lower) / se_bootstrap,
      ci_upper_critical = (ci_upper - original_fit$att) / se_bootstrap,
      t_stat_critical = scale_diagnostics$t_stat_critical,
      t_stat_critical_scale = scale_diagnostics$t_stat_critical_scale,
      ci_critical_scale = scale_diagnostics$ci_critical_scale,
      bootstrap_scale_t_stat = scale_diagnostics$bootstrap_scale_t_stat,
      pvalue = mean(abs(t_stats_valid) >= abs(t_stat_original)),
      rejection_rate = mean(abs(t_stats_valid) >= abs(t_stat_original)),
      n_clusters = as.integer(g),
      requested_n_bootstrap = weight_info$requested_n_bootstrap,
      actual_n_bootstrap = weight_info$actual_n_bootstrap,
      n_bootstrap = weight_info$actual_n_bootstrap,
      weight_type = as.character(weight_type),
      t_stat_original = t_stat_original,
      att_bootstrap = att_bootstrap,
      t_stats_bootstrap = t_stats_bootstrap,
      method = "native",
      impose_null = isTRUE(impose_null),
      full_enumeration = weight_info$full_enumeration,
      ci_method = ci_method,
      restricted_model = restricted_model,
      use_fwildclusterboot = FALSE
    ),
    class = "lwdid_wcb_result"
  )
}


.wcb_native_with_controls <- function(
    data, y_transformed, d, cluster_var, controls, n_bootstrap = 999L,
    weight_type = "rademacher", alpha = 0.05, impose_null = TRUE,
    full_enumeration = NULL, restricted_model = "with_controls"
) {
  y <- as.numeric(data[[y_transformed]])
  d_vec <- as.numeric(data[[d]])
  controls_info <- .wcb_validate_controls_matrix(data, controls)
  control_matrix <- controls_info$matrix
  cluster_index <- as.integer(as.factor(data[[cluster_var]]))
  g <- length(unique(cluster_index))

  design_matrix <- cbind("(Intercept)" = 1, D = d_vec, control_matrix)
  xtx <- crossprod(design_matrix)
  condition_number <- kappa(xtx, exact = TRUE)
  if (is.finite(condition_number) && condition_number > 1e10) {
    warning(
      lwdid_numerical_warning(
        detail = paste(
          "Design matrix condition number is large.",
          "Numerical accuracy may be reduced."
        ),
        source_function = ".wcb_native",
        condition_number = condition_number,
        call = NULL
      )
    )
  }

  if (!.wcb_treatment_is_identified(design_matrix, d_col_idx = 2L)) {
    return(
      .degenerate_wcb_result(
        att = NaN,
        original_se = NaN,
        requested_n_bootstrap = n_bootstrap,
        actual_n_bootstrap = n_bootstrap,
        weight_type = weight_type,
        impose_null = isTRUE(impose_null),
        full_enumeration = FALSE,
        ci_method = if (isTRUE(impose_null)) "percentile_t" else "percentile",
        n_clusters = as.integer(g),
        restricted_model = restricted_model
      )
    )
  }

  xtx_inv <- .wcb_solve_system(
    a = xtx,
    unavailable_message = "Native WCB controls design matrix is singular and MASS package unavailable."
  )
  original_fit <- .fit_wcb_linear_model(
    y = y,
    X = design_matrix,
    cluster_index = cluster_index,
    g = g,
    xtx_inv = xtx_inv,
    d_col_idx = 2L
  )
  if (.wcb_is_degenerate_standard_error(original_fit$se)) {
    return(
      .degenerate_wcb_result(
        att = original_fit$att,
        original_se = original_fit$se,
        requested_n_bootstrap = n_bootstrap,
        actual_n_bootstrap = n_bootstrap,
        weight_type = weight_type,
        impose_null = isTRUE(impose_null),
        full_enumeration = FALSE,
        ci_method = if (isTRUE(impose_null)) "percentile_t" else "percentile",
        n_clusters = as.integer(g),
        restricted_model = restricted_model
      )
    )
  }

  if (isTRUE(impose_null)) {
    if (identical(restricted_model, "intercept_only")) {
      restricted_matrix <- matrix(
        1,
        nrow = nrow(data),
        ncol = 1L,
        dimnames = list(NULL, "(Intercept)")
      )
    } else {
      restricted_matrix <- cbind("(Intercept)" = 1, control_matrix)
    }
    beta_restricted <- drop(
      .wcb_solve_system(
        a = crossprod(restricted_matrix),
        b = crossprod(restricted_matrix, y),
        unavailable_message = "Restricted WCB controls design matrix is singular and MASS package unavailable."
      )
    )
    fitted_base <- drop(restricted_matrix %*% beta_restricted)
    resid_base <- y - fitted_base
  } else {
    beta_restricted <- NULL
    fitted_base <- original_fit$fitted
    resid_base <- original_fit$residuals
  }

  .wcb_finalize_native_result(
    original_fit = original_fit,
    design_matrix = design_matrix,
    cluster_index = cluster_index,
    g = g,
    xtx_inv = xtx_inv,
    fitted_base = fitted_base,
    resid_base = resid_base,
    requested_n_bootstrap = n_bootstrap,
    weight_type = weight_type,
    alpha = alpha,
    impose_null = impose_null,
    full_enumeration = full_enumeration,
    restricted_model = restricted_model,
    restricted_coefficients = beta_restricted
  )
}


.fast_cluster_se <- function(
    residuals, X, obs_cluster_idx, G, XtX_inv, param_idx = 2L
) {
  residuals <- as.numeric(residuals)
  X <- as.matrix(X)
  obs_cluster_idx <- as.integer(obs_cluster_idx)
  XtX_inv <- as.matrix(XtX_inv)
  param_idx <- as.integer(param_idx)
  G <- as.integer(G)

  if (length(residuals) != nrow(X)) {
    stop("`residuals` must have length `nrow(X)`.", call. = FALSE)
  }
  if (length(obs_cluster_idx) != nrow(X)) {
    stop("`obs_cluster_idx` must have length `nrow(X)`.", call. = FALSE)
  }
  if (!all(dim(XtX_inv) == c(ncol(X), ncol(X)))) {
    stop("`XtX_inv` must be a square matrix matching `ncol(X)`.", call. = FALSE)
  }
  if (param_idx < 1L || param_idx > ncol(X)) {
    stop("`param_idx` is out of bounds for `X`.", call. = FALSE)
  }
  if (G <= 1L || nrow(X) <= ncol(X)) {
    return(NaN)
  }

  meat <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  for (cluster_rows in split(seq_len(nrow(X)), obs_cluster_idx)) {
    score_g <- drop(crossprod(X[cluster_rows, , drop = FALSE], residuals[cluster_rows]))
    meat <- meat + tcrossprod(score_g)
  }

  correction <- (G / (G - 1)) * ((nrow(X) - 1) / (nrow(X) - ncol(X)))
  var_beta <- correction * XtX_inv %*% meat %*% XtX_inv

  sqrt(max(as.numeric(var_beta[param_idx, param_idx]), 0))
}


.generate_weights <- function(n_clusters, weight_type = "rademacher") {
  n_clusters <- as.integer(n_clusters)
  if (n_clusters <= 0L) {
    stop("`n_clusters` must be positive.", call. = FALSE)
  }

  if (identical(weight_type, "rademacher")) {
    return(sample(c(-1, 1), size = n_clusters, replace = TRUE))
  }

  if (identical(weight_type, "mammen")) {
    sqrt5 <- sqrt(5)
    phi <- (sqrt5 + 1) / 2
    p_small <- (sqrt5 + 1) / (2 * sqrt5)
    draws <- stats::runif(n_clusters)

    return(ifelse(draws < p_small, -(phi - 1), phi))
  }

  if (identical(weight_type, "webb")) {
    webb_support <- c(-sqrt(3 / 2), -1, -sqrt(1 / 2), sqrt(1 / 2), 1, sqrt(3 / 2))

    return(sample(webb_support, size = n_clusters, replace = TRUE))
  }

  stop(
    sprintf(
      "Unsupported weight_type '%s'. Must be one of: 'rademacher', 'mammen', 'webb'.",
      weight_type
    ),
    call. = FALSE
  )
}


.generate_all_rademacher <- function(n_clusters) {
  n_clusters <- as.integer(n_clusters)
  if (n_clusters <= 0L) {
    stop("`n_clusters` must be positive.", call. = FALSE)
  }

  weights <- expand.grid(
    rep(list(c(-1, 1)), n_clusters),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  data.matrix(weights)
}


.degenerate_wcb_result <- function(
    att, original_se = NaN, requested_n_bootstrap, actual_n_bootstrap,
    weight_type, impose_null, full_enumeration, ci_method, n_clusters,
    restricted_model = "with_controls"
) {
  original_se_value <- unname(as.numeric(original_se))
  if (!is.finite(original_se_value) || original_se_value <= 0) {
    original_se_value <- NaN
  }

  structure(
    list(
      att = unname(as.numeric(att)),
      original_se = original_se_value,
      se_bootstrap = NaN,
      ci_lower = NaN,
      ci_upper = NaN,
      ci_lower_critical = NaN,
      ci_upper_critical = NaN,
      t_stat_critical = NaN,
      t_stat_critical_scale = "original_t",
      ci_critical_scale = "bootstrap_se",
      bootstrap_scale_t_stat = NaN,
      pvalue = NaN,
      rejection_rate = NaN,
      n_clusters = as.integer(n_clusters),
      requested_n_bootstrap = as.integer(requested_n_bootstrap),
      actual_n_bootstrap = as.integer(actual_n_bootstrap),
      n_bootstrap = as.integer(actual_n_bootstrap),
      weight_type = as.character(weight_type),
      t_stat_original = NaN,
      att_bootstrap = numeric(0),
      t_stats_bootstrap = numeric(0),
      method = "native",
      impose_null = isTRUE(impose_null),
      full_enumeration = isTRUE(full_enumeration),
      ci_method = as.character(ci_method),
      restricted_model = restricted_model,
      use_fwildclusterboot = FALSE
    ),
    class = "lwdid_wcb_result"
  )
}


.wcb_is_degenerate_standard_error <- function(
    se,
    treat_near_zero_as_degenerate = FALSE
) {
  if (!is.finite(se)) {
    return(TRUE)
  }

  threshold <- if (isTRUE(treat_near_zero_as_degenerate)) {
    LWDID_VARIANCE_THRESHOLD
  } else {
    0
  }

  as.numeric(se) <= threshold
}


.wcb_scale_diagnostics <- function(t_stats_valid, att, se_bootstrap, alpha = 0.05) {
  valid_t <- abs(as.numeric(t_stats_valid))
  valid_t <- valid_t[is.finite(valid_t)]

  list(
    t_stat_critical = if (length(valid_t) == 0L) {
      NaN
    } else {
      as.numeric(
        stats::quantile(valid_t, probs = 1 - alpha, names = FALSE, type = 7)
      )
    },
    t_stat_critical_scale = "original_t",
    ci_critical_scale = "bootstrap_se",
    bootstrap_scale_t_stat = if (is.finite(se_bootstrap) && se_bootstrap > 0) {
      abs(as.numeric(att) / as.numeric(se_bootstrap))
    } else {
      NaN
    }
  )
}


.compute_bootstrap_pvalue_at_null <- function(
    data, y_transformed, d, cluster_var, controls = NULL, null_value,
    att_original, se_original, n_bootstrap, weight_type, seed = NULL,
    full_enumeration = NULL
) {
  if (.wcb_is_degenerate_standard_error(
    se_original,
    treat_near_zero_as_degenerate = TRUE
  )) {
    return(NaN)
  }

  data_adj <- data
  data_adj[[y_transformed]] <- data[[y_transformed]] - null_value * data[[d]]
  t_orig <- (att_original - null_value) / se_original

  result <- wild_cluster_bootstrap(
    data = data_adj,
    y_transformed = y_transformed,
    d = d,
    cluster_var = cluster_var,
    controls = controls,
    n_bootstrap = n_bootstrap,
    weight_type = weight_type,
    seed = seed,
    impose_null = TRUE,
    full_enumeration = full_enumeration,
    restricted_model = "intercept_only",
    use_fwildclusterboot = FALSE
  )

  valid_t <- result$t_stats_bootstrap[is.finite(result$t_stats_bootstrap)]
  if (length(valid_t) == 0L) {
    return(NaN)
  }

  mean(abs(valid_t) >= abs(t_orig))
}


#' Wild cluster bootstrap test-inversion confidence interval
#'
#' This public native slice inverts bootstrap hypothesis tests on the native
#' path. It currently uses the Python-compatible intercept-only restricted
#' model, while `fwildclusterboot` integration remains deferred.
#'
#' @inheritParams wild_cluster_bootstrap
#' @param grid_points Integer scalar giving the coarse search grid size.
#' @param ci_tol Numeric scalar giving the root-finding tolerance.
#'
#' @return An object of class `lwdid_wcb_result`.
#' @export
wild_cluster_bootstrap_test_inversion <- function(
    data, y_transformed, d, cluster_var, controls = NULL,
    n_bootstrap = 999L, weight_type = "rademacher", alpha = 0.05,
    seed = NULL, grid_points = 25L, ci_tol = 0.01
) {
  original_result <- wild_cluster_bootstrap(
    data = data,
    y_transformed = y_transformed,
    d = d,
    cluster_var = cluster_var,
    controls = controls,
    n_bootstrap = n_bootstrap,
    weight_type = weight_type,
    alpha = alpha,
    seed = seed,
    impose_null = TRUE,
    full_enumeration = FALSE,
    restricted_model = "intercept_only",
    use_fwildclusterboot = FALSE
  )

  att_original <- original_result$att
  se_original <- original_result$original_se
  if (.wcb_is_degenerate_standard_error(
    se_original,
    treat_near_zero_as_degenerate = TRUE
  )) {
    return(
      .degenerate_wcb_result(
        att = att_original,
        original_se = se_original,
        requested_n_bootstrap = as.integer(n_bootstrap),
        actual_n_bootstrap = as.integer(n_bootstrap),
        weight_type = weight_type,
        impose_null = TRUE,
        full_enumeration = FALSE,
        ci_method = "test_inversion",
        n_clusters = as.integer(original_result$n_clusters),
        restricted_model = "intercept_only"
      )
    )
  }

  pvalue_at_theta <- function(theta) {
    .compute_bootstrap_pvalue_at_null(
      data = data,
      y_transformed = y_transformed,
      d = d,
      cluster_var = cluster_var,
      controls = controls,
      null_value = theta,
      att_original = att_original,
      se_original = se_original,
      n_bootstrap = n_bootstrap,
      weight_type = weight_type,
      seed = seed,
      full_enumeration = FALSE
    )
  }

  search_range <- 4 * se_original
  theta_grid <- seq(
    att_original - search_range,
    att_original + search_range,
    length.out = as.integer(grid_points)
  )
  pvalues_grid <- vapply(theta_grid, pvalue_at_theta, numeric(1))
  ci_mask <- pvalues_grid >= alpha

  if (!any(ci_mask)) {
    ci_lower <- att_original - 1.96 * se_original
    ci_upper <- att_original + 1.96 * se_original
  } else {
    ci_lower_approx <- min(theta_grid[ci_mask])
    ci_upper_approx <- max(theta_grid[ci_mask])
    pvalue_minus_alpha <- function(theta) {
      pvalue_at_theta(theta) - alpha
    }

    ci_lower <- tryCatch(
      stats::uniroot(
        pvalue_minus_alpha,
        interval = c(att_original - search_range, att_original),
        tol = ci_tol
      )$root,
      error = function(...) ci_lower_approx
    )
    ci_upper <- tryCatch(
      stats::uniroot(
        pvalue_minus_alpha,
        interval = c(att_original, att_original + search_range),
        tol = ci_tol
      )$root,
      error = function(...) ci_upper_approx
    )
  }

  se_bootstrap <- (ci_upper - ci_lower) / (2 * 1.96)
  structure(
    list(
      att = att_original,
      original_se = se_original,
      se_bootstrap = se_bootstrap,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_lower_critical = (att_original - ci_lower) / se_bootstrap,
      ci_upper_critical = (ci_upper - att_original) / se_bootstrap,
      t_stat_critical = NaN,
      t_stat_critical_scale = "original_t",
      ci_critical_scale = "bootstrap_se",
      bootstrap_scale_t_stat = abs(att_original / se_bootstrap),
      pvalue = pvalue_at_theta(0),
      rejection_rate = pvalue_at_theta(0),
      n_clusters = as.integer(original_result$n_clusters),
      requested_n_bootstrap = as.integer(n_bootstrap),
      actual_n_bootstrap = as.integer(n_bootstrap),
      n_bootstrap = as.integer(n_bootstrap),
      weight_type = as.character(weight_type),
      t_stat_original = att_original / se_original,
      att_bootstrap = numeric(0),
      t_stats_bootstrap = numeric(0),
      method = "native",
      impose_null = TRUE,
      full_enumeration = FALSE,
      ci_method = "test_inversion",
      restricted_model = "intercept_only",
      use_fwildclusterboot = FALSE
    ),
    class = "lwdid_wcb_result"
  )
}


#' Print a wild cluster bootstrap result
#'
#' @param x An object of class `lwdid_wcb_result`.
#' @param digits Number of printed digits.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.lwdid_wcb_result <- function(x, digits = 4L, ...) {
  cat("Wild Cluster Bootstrap Result\n")
  cat(sprintf("ATT: %s\n", formatC(x$att, format = "f", digits = digits)))
  cat(
    sprintf(
      "P-value: %s | CI: [%s, %s]\n",
      formatC(x$pvalue, format = "f", digits = digits),
      formatC(x$ci_lower, format = "f", digits = digits),
      formatC(x$ci_upper, format = "f", digits = digits)
    )
  )
  cat(
    sprintf(
      "Clusters: %d | Requested bootstrap: %d | Actual bootstrap: %d\n",
      x$n_clusters,
      x$requested_n_bootstrap,
      x$actual_n_bootstrap
    )
  )
  cat(
    sprintf(
      "Weight type: %s | Full enumeration: %s\n",
      x$weight_type,
      if (isTRUE(x$full_enumeration)) "TRUE" else "FALSE"
    )
  )

  invisible(x)
}


#' Summarise a wild cluster bootstrap result
#'
#' @param object An object of class `lwdid_wcb_result`.
#' @param ... Unused.
#'
#' @return `object`, invisibly.
#' @export
summary.lwdid_wcb_result <- function(object, ...) {
  print(object, ...)
  invisible(object)
}
