library(testthat)

resolve_package_root <- function() {
  candidates <- character(0)

  namespace_root <- tryCatch(
    getNamespaceInfo(asNamespace("lwdid"), "path"),
    error = function(...) NULL
  )
  if (!is.null(namespace_root) && nzchar(namespace_root)) {
    candidates <- c(candidates, namespace_root)
  }

  test_path_candidate <- tryCatch(
    testthat::test_path("..", ".."),
    error = function(...) NULL
  )
  if (!is.null(test_path_candidate) && nzchar(test_path_candidate)) {
    candidates <- c(candidates, test_path_candidate)
  }

  candidates <- c(candidates, ".")
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))

  existing <- candidates[file.exists(file.path(candidates, "NAMESPACE"))]
  if (length(existing) > 0L) {
    return(normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE))
  }

  normalizePath(candidates[[1L]], winslash = "/", mustWork = FALSE)
}

read_generated_namespace <- function() {
  pkg_root <- resolve_package_root()
  parseNamespaceFile(basename(pkg_root), dirname(pkg_root))
}

has_registered_s3_method <- function(ns_info, generic, class_name) {
  s3_methods <- ns_info$S3methods
  if (is.null(s3_methods) || length(s3_methods) == 0L) {
    return(FALSE)
  }

  any(s3_methods[, 1] == generic & s3_methods[, 2] == class_name)
}

test_that("E8-03 namespace contract exports lwdid_sensitivity", {
  ns_info <- read_generated_namespace()

  expect_true("lwdid_sensitivity" %in% ns_info$exports)
})

test_that("E8-03 namespace contract registers comprehensive print and summary methods", {
  ns_info <- read_generated_namespace()

  expect_true(has_registered_s3_method(
    ns_info,
    "print",
    "lwdid_sensitivity_comprehensive"
  ))
  expect_true(has_registered_s3_method(
    ns_info,
    "summary",
    "lwdid_sensitivity_comprehensive"
  ))
})
