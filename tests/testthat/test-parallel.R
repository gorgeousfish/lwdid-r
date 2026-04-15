# test-parallel.R — Tests for parallel computation infrastructure
# Epic 13: Performance/Parallel

# =========================================================================
# safe_iterate tests
# =========================================================================

test_that("safe_iterate: success returns structured result", {
  result <- lwdid:::safe_iterate(function() 42)
  expect_equal(result$status, "success")
  expect_equal(result$result, 42)
  expect_null(result$error)
  expect_length(result$warnings, 0L)
})

test_that("safe_iterate: error returns failed status", {
  result <- lwdid:::safe_iterate(function() stop("test error"))
  expect_equal(result$status, "failed")
  expect_null(result$result)
  expect_match(result$error, "test error")
})

test_that("safe_iterate: warnings captured without stopping", {
  result <- lwdid:::safe_iterate(function() {
    warning("warn1")
    warning("warn2")
    99
  })
  expect_equal(result$status, "success")
  expect_equal(result$result, 99)
  expect_length(result$warnings, 2L)
  expect_match(result$warnings[1], "warn1")
})

# =========================================================================
# run_parallel tests (sequential mode)
# =========================================================================

test_that("run_parallel: sequential mode returns correct results", {
  res <- run_parallel(1:5, function(x) x^2, parallel = FALSE)
  expect_equal(unlist(res), c(1, 4, 9, 16, 25))
})

test_that("run_parallel: empty input returns empty list", {
  res <- run_parallel(integer(0), function(x) x, parallel = FALSE)
  expect_length(res, 0L)
  expect_true(is.list(res))
})

test_that("run_parallel: single element works", {
  res <- run_parallel(1L, function(x) x^2, parallel = FALSE)
  expect_length(res, 1L)
  expect_equal(res[[1]], 1L)
})

test_that("run_parallel: extra args passed via ...", {
  multiply_fn <- function(b, factor) b * factor
  res <- run_parallel(1:3, multiply_fn, factor = 10, parallel = FALSE)
  expect_equal(unlist(res), c(10, 20, 30))
})

test_that("run_parallel: multiple extra args", {
  compute_fn <- function(b, offset, scale) (b + offset) * scale
  res <- run_parallel(1:3, compute_fn, offset = 100, scale = 2,
                      parallel = FALSE)
  expect_equal(unlist(res), c(202, 204, 206))
})

# =========================================================================
# run_parallel fault tolerance tests
# =========================================================================

test_that("run_parallel: single failure does not stop loop", {
  res <- run_parallel(
    1:10,
    function(b) {
      if (b == 5) stop("deliberate error")
      b^2
    },
    parallel = FALSE,
    fail_threshold = 1.0
  )
  expect_length(res, 9L)
  expect_equal(unlist(res), (1:10)[-5]^2)
})

test_that("run_parallel: >50% failure throws error", {
  expect_error(
    run_parallel(
      1:10,
      function(b) {
        if (b <= 6) stop("fail")
        b
      },
      parallel = FALSE,
      fail_threshold = 1.0
    ),
    "exceeds 50%"
  )
})

test_that("run_parallel: failure above threshold triggers warning", {
  expect_warning(
    run_parallel(
      1:10,
      function(b) {
        if (b <= 2) stop("fail")
        b
      },
      parallel = FALSE,
      fail_threshold = 0.1
    ),
    "iterations failed"
  )
})

# =========================================================================
# setup_parallel user-facing tests
# =========================================================================

test_that("setup_parallel: sequential strategy", {
  skip_if_not_installed("future")
  old <- future::plan("list")
  on.exit(future::plan(old), add = TRUE)

  expect_message(setup_parallel(strategy = "sequential"),
                 "sequential")
})

test_that("setup_parallel: multisession strategy", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()
  old <- future::plan("list")
  on.exit(future::plan(old), add = TRUE)

  expect_message(setup_parallel(n_workers = 2, strategy = "multisession"),
                 "parallel computation enabled")
})

# =========================================================================
# .setup_parallel_internal tests
# =========================================================================

test_that(".setup_parallel_internal: parallel=FALSE returns NULL", {
  res <- lwdid:::.setup_parallel_internal(parallel = FALSE)
  expect_null(res)
})

test_that(".setup_parallel_internal: n_cores < 1 gives error", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  expect_error(
    lwdid:::.setup_parallel_internal(parallel = TRUE, n_cores = 0L),
    "n_cores must be >= 1"
  )
})

test_that(".setup_parallel_internal: n_cores=1 returns NULL (no plan change)", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  res <- lwdid:::.setup_parallel_internal(parallel = TRUE, n_cores = 1L)
  expect_null(res)
})

test_that(".setup_parallel_internal: returns cleanup function", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()
  old <- future::plan("list")
  on.exit(future::plan(old), add = TRUE)

  future::plan(future::sequential)
  cleanup <- lwdid:::.setup_parallel_internal(parallel = TRUE, n_cores = 2L)
  expect_true(is.function(cleanup))
  # Plan should now be non-sequential
  expect_false(inherits(future::plan(), "sequential"))
  # Cleanup restores
  cleanup()
  expect_true(inherits(future::plan(), "sequential"))
})

test_that(".setup_parallel_internal: respects existing non-sequential plan", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()
  old <- future::plan("list")
  on.exit(future::plan(old), add = TRUE)

  future::plan(future::multisession, workers = 2L)
  # Should return NULL (no change) since plan already non-sequential
  res <- lwdid:::.setup_parallel_internal(parallel = TRUE, n_cores = 4L)
  expect_null(res)
})

# =========================================================================
# run_parallel with actual parallel execution
# =========================================================================

test_that("run_parallel: parallel=TRUE works with future", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()

  res <- run_parallel(
    1:10, function(x) x^2,
    parallel = TRUE, n_cores = 2L,
    future.seed = NULL
  )
  expect_equal(sort(unlist(res)), (1:10)^2)
})

test_that("run_parallel: deterministic with future.seed=TRUE", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()

  set.seed(123)
  res1 <- run_parallel(
    1:20, function(b) rnorm(1),
    parallel = TRUE, n_cores = 2L,
    future.seed = TRUE
  )
  set.seed(123)
  res2 <- run_parallel(
    1:20, function(b) rnorm(1),
    parallel = TRUE, n_cores = 2L,
    future.seed = TRUE
  )
  expect_equal(unlist(res1), unlist(res2))
})

# =========================================================================
# Staggered integration test
# =========================================================================

test_that("staggered estimation works with run_parallel refactor", {
  data(castle, package = "lwdid")
  dt <- data.table::as.data.table(castle)

  result <- suppressWarnings(lwdid(
    data = dt,
    y = "lhomicide",
    gvar = "effyear",
    ivar = "sid",
    tvar = "year",
    rolling = "demean",
    aggregate = "overall"
  ))
  expect_s3_class(result, "lwdid_result")
  expect_true(is.numeric(result$att))
  expect_true(is.finite(result$att))
})
