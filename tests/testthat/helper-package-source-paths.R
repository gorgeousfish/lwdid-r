resolve_package_source_root <- function() {
  candidates <- character(0)

  source_copy_candidate <- tryCatch(
    testthat::test_path("..", "..", "00_pkg_src", "lwdid"),
    error = function(...) NULL
  )
  if (!is.null(source_copy_candidate) && nzchar(source_copy_candidate)) {
    candidates <- c(candidates, source_copy_candidate)
  }

  candidates <- c(candidates, "/Users/cxy/Desktop/lwdid_r/lwdid-r")

  test_path_candidate <- tryCatch(
    testthat::test_path("..", ".."),
    error = function(...) NULL
  )
  if (!is.null(test_path_candidate) && nzchar(test_path_candidate)) {
    candidates <- c(candidates, test_path_candidate)
  }

  namespace_root <- tryCatch(
    getNamespaceInfo(asNamespace("lwdid"), "path"),
    error = function(...) NULL
  )
  if (!is.null(namespace_root) && nzchar(namespace_root)) {
    candidates <- c(candidates, namespace_root)
  }

  candidates <- c(candidates, ".")
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))

  existing <- candidates[file.exists(file.path(candidates, "NAMESPACE"))]
  if (length(existing) > 0L) {
    return(normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE))
  }

  normalizePath(candidates[[1L]], winslash = "/", mustWork = FALSE)
}

resolve_package_source_file <- function(...) {
  file.path(resolve_package_source_root(), ...)
}
