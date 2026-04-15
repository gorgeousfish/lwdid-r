resolve_parity_fixture_path <- function(filename) {
  candidates <- character(0)
  repo_root_env <- Sys.getenv("LWDID_REPO_ROOT", unset = "")

  repo_path_candidate <- tryCatch(
    testthat::test_path(
      "..", "..", "..",
      "_automation", "test-artifacts", "parity", filename
    ),
    error = function(...) NULL
  )
  if (!is.null(repo_path_candidate) && nzchar(repo_path_candidate)) {
    candidates <- c(candidates, repo_path_candidate)
  }

  test_path_candidate <- tryCatch(
    testthat::test_path("..", "p", filename),
    error = function(...) NULL
  )
  if (!is.null(test_path_candidate) && nzchar(test_path_candidate)) {
    candidates <- c(candidates, test_path_candidate)
  }

  if (nzchar(repo_root_env)) {
    candidates <- c(
      candidates,
      file.path(repo_root_env, "_automation", "test-artifacts", "parity", filename)
    )
  }

  candidates <- c(
    candidates,
    file.path(
      "/Users/cxy/Desktop/lwdid_r",
      "_automation", "test-artifacts", "parity", filename
    )
  )

  namespace_root <- tryCatch(
    getNamespaceInfo(asNamespace("lwdid"), "path"),
    error = function(...) NULL
  )
  if (!is.null(namespace_root) && nzchar(namespace_root)) {
    candidates <- c(
      candidates,
      file.path(namespace_root, "tests", "p", filename),
      file.path(
        normalizePath(file.path(namespace_root, ".."), mustWork = FALSE),
        "_automation", "test-artifacts", "parity", filename
      )
    )
  }

  candidates <- c(
    candidates,
    file.path("..", "_automation", "test-artifacts", "parity", filename)
  )
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))

  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) {
    return(existing[[1L]])
  }

  candidates[[1L]]
}

skip_if_no_parity_fixture <- function(filename) {
  path <- resolve_parity_fixture_path(filename)
  if (!file.exists(path)) {
    testthat::skip(paste("Parity fixture not available:", filename))
  }
  invisible(path)
}

parity_fixtures_available <- function() {
  test_p <- tryCatch(
    testthat::test_path("..", "p"),
    error = function(...) NULL
  )
  if (!is.null(test_p) && dir.exists(test_p)) return(TRUE)
  repo_p <- tryCatch(
    testthat::test_path("..", "..", "..", "_automation", "test-artifacts", "parity"),
    error = function(...) NULL
  )
  if (!is.null(repo_p) && dir.exists(repo_p)) return(TRUE)
  dir.exists("/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity")
}

skip_if_no_parity_fixtures <- function() {
  if (!parity_fixtures_available()) {
    testthat::skip("Parity fixtures not available (tests/p excluded from build)")
  }
}
