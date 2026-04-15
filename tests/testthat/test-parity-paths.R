library(testthat)

test_that("parity fixture resolver prefers live repo assets over tests/p snapshots", {
  filename <- "20260325-qa-parity-e8-06-trend-public-staggered-fallback.json"
  repo_path <- testthat::test_path(
    "..", "..", "..",
    "_automation", "test-artifacts", "parity", filename
  )
  stale_snapshot_path <- testthat::test_path("..", "p", filename)

  skip_if_not(file.exists(repo_path), paste("missing repo parity oracle:", repo_path))
  skip_if(
    file.exists(stale_snapshot_path),
    paste("refusing to overwrite existing snapshot:", stale_snapshot_path)
  )

  dir.create(dirname(stale_snapshot_path), recursive = TRUE, showWarnings = FALSE)
  writeLines('{"stale_snapshot": true}', stale_snapshot_path)
  on.exit(unlink(stale_snapshot_path), add = TRUE)

  expect_identical(
    normalizePath(resolve_parity_fixture_path(filename), mustWork = FALSE),
    normalizePath(repo_path, mustWork = FALSE)
  )
})
