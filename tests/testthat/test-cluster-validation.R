# ============================================================================
# Tests for Cluster Pre-Validation and Nesting (Task E3-05.6)
# Story E3-05: VCE Integration — Cluster Validation Tests
#
# Tests cluster validation logic in .estimate_common_timing() (Step 1c)
# and graded warnings in compute_cluster_vce().
#
# Test groups:
#   CV-01: cluster_var column not found
#   CV-02: cluster_var column contains NA
#   CV-03: Cluster nesting violation
#   CV-04: Correct nesting passes
#   CV-05: Graded warning G<10
#   CV-06: Graded warning 10<=G<20
#   CV-07: G=10 boundary
#   CV-08: G>=20 no small-sample warning
#   CV-09: Cluster size imbalance CV>1.0
#   CV-10: Cluster size balanced CV<=1.0
#   CV-11: Character cluster variable
#   CV-12: Numeric cluster variable (same VCE as character)
# ============================================================================

library(testthat)
library(data.table)

# ============================================================================
# Group 1: Pre-validation in .estimate_common_timing() via lwdid()
# CV-01 through CV-04 test the cluster validation in Step 1c of
# .estimate_common_timing(), which runs on the original panel data
# before transformation.
# ============================================================================

# CV-01: cluster_var column not found → lwdid_error
# validate_inputs() catches nonexistent columns first with
# lwdid_missing_column (a subclass of lwdid_error).
test_that("CV-01: cluster_var column not found triggers error", {
  data(smoking)
  expect_error(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "nonexistent_column"),
    class = "lwdid_error"
  )
  # Verify the error message mentions the nonexistent column
  err <- tryCatch(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "nonexistent_column"),
    error = function(e) e
  )
  expect_true(grepl("nonexistent_column", conditionMessage(err)))
  expect_true(grepl("not found", conditionMessage(err)))
})

# CV-02: cluster_var column contains NA → lwdid_invalid_parameter
test_that("CV-02: cluster_var with NAs triggers lwdid_invalid_parameter", {
  data(smoking)
  dt <- data.table::copy(data.table::as.data.table(smoking))
  dt[, cluster_col := state]
  dt[1:3, cluster_col := NA]
  expect_error(
    suppressWarnings(
      lwdid(dt, y = "lcigsale", ivar = "state", tvar = "year",
            d = "d", post = "post", rolling = "demean",
            vce = "cluster", cluster_var = "cluster_col")
    ),
    class = "lwdid_invalid_parameter"
  )
  # Verify error message mentions NA count
  err <- tryCatch(
    suppressWarnings(
      lwdid(dt, y = "lcigsale", ivar = "state", tvar = "year",
            d = "d", post = "post", rolling = "demean",
            vce = "cluster", cluster_var = "cluster_col")
    ),
    error = function(e) e
  )
  expect_true(grepl("3 NA", conditionMessage(err)))
})

# CV-03: Cluster nesting violation → lwdid_invalid_parameter
# When a unit belongs to multiple clusters across time periods.
test_that("CV-03: cluster nesting violation triggers lwdid_invalid_parameter", {
  set.seed(42)
  n_units <- 10L
  t_periods <- 5L
  dt <- data.table::CJ(unit = 1:n_units, year = 1:t_periods)
  dt[, d := as.integer(unit <= 3)]
  dt[, post := as.integer(year >= 4)]
  dt[, y := rnorm(.N)]
  # Assign clusters: units 1-5 -> "A", units 6-10 -> "B"
  dt[, cluster := ifelse(unit <= 5, "A", "B")]
  # Make unit 1 violate nesting: cluster=A in years 1-3, cluster=B in 4-5
  dt[unit == 1 & year >= 4, cluster := "B"]

  expect_error(
    suppressWarnings(
      lwdid(dt, y = "y", ivar = "unit", tvar = "year",
            d = "d", post = "post", rolling = "demean",
            vce = "cluster", cluster_var = "cluster")
    ),
    class = "lwdid_invalid_parameter"
  )

  # Verify error message mentions nesting violation and unit examples
  err <- tryCatch(
    suppressWarnings(
      lwdid(dt, y = "y", ivar = "unit", tvar = "year",
            d = "d", post = "post", rolling = "demean",
            vce = "cluster", cluster_var = "cluster")
    ),
    error = function(e) e
  )
  expect_true(grepl("nesting", conditionMessage(err), ignore.case = TRUE))
})

# CV-04: Correct nesting → normal computation
test_that("CV-04: correct nesting produces valid results", {
  data(smoking)
  result <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "state")
  )
  # Should succeed without error
  expect_true(!is.null(result$att))
  expect_true(!is.null(result$se_att))
  expect_true(result$se_att > 0)
  # ATT should be finite and numerically reasonable
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se_att))
})

# ============================================================================
# Group 2: Graded cluster warnings from compute_cluster_vce()
# CV-05 through CV-08 test the graded warning thresholds.
# These test estimate_ra_common() directly since the warnings come
# from compute_cluster_vce() which is called inside it.
# ============================================================================

# CV-05: G<10 triggers strong lwdid_small_sample warning
test_that("CV-05: G<10 triggers strong small-sample warning", {
  set.seed(55)
  g <- 5L
  n_per <- 6L
  n <- g * n_per
  cluster <- rep(seq_len(g), each = n_per)
  d <- c(rep(1L, 10), rep(0L, 20))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # Check for strong small-sample warning
  small_sample_warns <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_caught
  )
  expect_true(length(small_sample_warns) >= 1L)
  # Message should indicate "highly unreliable"
  msgs <- vapply(small_sample_warns, conditionMessage, character(1))
  expect_true(any(grepl("highly unreliable", msgs, ignore.case = TRUE)))
  # Should mention G=5
  expect_true(any(grepl("G=5", msgs)))
})

# CV-06: 10<=G<20 triggers informational lwdid_small_sample warning
test_that("CV-06: 10<=G<20 triggers informational warning", {
  set.seed(66)
  g <- 15L
  n_per <- 4L
  n <- g * n_per
  cluster <- rep(seq_len(g), each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  small_sample_warns <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_caught
  )
  expect_true(length(small_sample_warns) >= 1L)
  msgs <- vapply(small_sample_warns, conditionMessage, character(1))
  # Should NOT contain "highly unreliable"
  expect_false(any(grepl("highly unreliable", msgs, ignore.case = TRUE)))
  # Should suggest Wild Cluster Bootstrap
  expect_true(any(grepl("Wild Cluster Bootstrap|WCB|bootstrap",
                        msgs, ignore.case = TRUE)))
  # Should mention G=15
  expect_true(any(grepl("G=15", msgs)))
})

# CV-07: G=10 boundary triggers informational warning (NOT strong)
test_that("CV-07: G=10 triggers informational warning, not strong", {
  set.seed(77)
  g <- 10L
  n_per <- 4L
  n <- g * n_per
  cluster <- rep(seq_len(g), each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  small_sample_warns <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_caught
  )
  expect_true(length(small_sample_warns) >= 1L)
  msgs <- vapply(small_sample_warns, conditionMessage, character(1))
  # G=10 is in [10, 20) range: informational, NOT strong
  expect_false(any(grepl("highly unreliable", msgs, ignore.case = TRUE)))
  # Should suggest Bootstrap
  expect_true(any(grepl("Bootstrap", msgs, ignore.case = TRUE)))
  # Should mention G=10
  expect_true(any(grepl("G=10", msgs)))
})

# CV-08: G>=20 no small-sample cluster warning
test_that("CV-08: G>=20 emits no small-sample cluster warning", {
  set.seed(88)
  g <- 20L
  n_per <- 4L
  n <- g * n_per
  cluster <- rep(seq_len(g), each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # Filter for cluster-count related warnings (not imbalance)
  cluster_count_warns <- Filter(
    function(w) {
      msg <- conditionMessage(w)
      inherits(w, "lwdid_small_sample") &&
      (grepl("few cluster", msg, ignore.case = TRUE) ||
       grepl("Very few", msg, ignore.case = TRUE) ||
       grepl("G=", msg)) &&
      !grepl("imbalance", msg, ignore.case = TRUE) &&
      !grepl("unbalanced", msg, ignore.case = TRUE) &&
      !grepl("CV=", msg)
    },
    warnings_caught
  )
  expect_equal(length(cluster_count_warns), 0L)
})

# ============================================================================
# Group 3: Cluster size imbalance warnings
# CV-09 and CV-10 test the CV > 1.0 imbalance check.
# Use G>=20 to avoid mixing with small-sample warnings.
# ============================================================================

# CV-09: Cluster size imbalance CV>1.0 triggers warning
test_that("CV-09: cluster size imbalance CV>1.0 triggers warning", {
  set.seed(99)
  # 25 clusters: 1 large (50 obs) + 24 small (2 obs each) = 98 obs
  # CV = sd(sizes)/mean(sizes) >> 1.0
  cluster <- c(rep(1L, 50), rep(2:25, each = 2))
  n <- length(cluster)
  d <- c(rep(1L, 40), rep(0L, n - 40))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  # Verify CV > 1.0
  sizes <- as.integer(table(cluster))
  cv <- sd(sizes) / mean(sizes)
  expect_true(cv > 1.0)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # Should have imbalance warning
  small_sample_warns <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_caught
  )
  expect_true(length(small_sample_warns) >= 1L)
  msgs <- vapply(small_sample_warns, conditionMessage, character(1))
  expect_true(any(grepl("imbalance|unequal|unbalanced",
                        msgs, ignore.case = TRUE)))
  # Should mention CV value
  expect_true(any(grepl("CV", msgs, ignore.case = TRUE)))
})

# CV-10: Cluster size balanced CV<=1.0 no imbalance warning
test_that("CV-10: balanced cluster sizes emit no imbalance warning", {
  set.seed(100)
  g <- 20L
  n_per <- 4L
  n <- g * n_per
  cluster <- rep(seq_len(g), each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  # Verify CV = 0 (perfectly balanced)
  sizes <- as.integer(table(cluster))
  cv <- sd(sizes) / mean(sizes)
  expect_equal(cv, 0)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # No imbalance warnings
  imbalance_warns <- Filter(
    function(w) {
      msg <- conditionMessage(w)
      grepl("imbalance", msg, ignore.case = TRUE) ||
      grepl("unbalanced", msg, ignore.case = TRUE) ||
      grepl("CV=", msg)
    },
    warnings_caught
  )
  expect_equal(length(imbalance_warns), 0L)
})

# ============================================================================
# Group 4: Character and numeric cluster variable types
# CV-11 and CV-12 verify both types work and produce identical VCE.
# ============================================================================

# CV-11: Character cluster variable works correctly
test_that("CV-11: character cluster variable works correctly", {
  set.seed(111)
  g <- 20L
  n_per <- 4L
  n <- g * n_per
  state_names <- paste0("S", sprintf("%02d", seq_len(g)))
  cluster_char <- rep(state_names, each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  res <- suppressWarnings(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster_char)
  )
  expect_true(is.finite(res$att))
  expect_true(is.finite(res$se))
  expect_equal(res$vce_type, "cluster")
  expect_equal(res$n_clusters, g)
})

# CV-12: Numeric and character cluster produce same VCE
test_that("CV-12: numeric and character cluster produce same VCE", {
  set.seed(112)
  g <- 20L
  n_per <- 4L
  n <- g * n_per
  cluster_num <- rep(seq_len(g), each = n_per)
  cluster_char <- rep(paste0("C", seq_len(g)), each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  res_num <- suppressWarnings(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster_num)
  )
  res_char <- suppressWarnings(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster_char)
  )

  # Both should produce valid results
  expect_true(is.finite(res_num$att))
  expect_true(is.finite(res_char$att))
  expect_equal(res_num$vce_type, "cluster")
  expect_equal(res_char$vce_type, "cluster")
  expect_equal(res_num$n_clusters, g)
  expect_equal(res_char$n_clusters, g)

  # ATT should be identical (same data, same model)
  expect_equal(res_num$att, res_char$att, tolerance = 1e-12)

  # VCE matrices should be identical
  # sandwich::vcovCL internally converts to factor, so ordering may
  # differ for character vs numeric. But since both use the same
  # underlying grouping, the VCE should match.
  expect_equal(res_num$vcov, res_char$vcov, tolerance = 1e-12)
  expect_equal(res_num$se, res_char$se, tolerance = 1e-12)
})

# ============================================================================
# Group 5: Boundary and edge case tests
# ============================================================================

# CV-05b: G=9 boundary (just below 10) triggers strong warning
test_that("CV-05b: G=9 triggers strong warning (G<10 path)", {
  set.seed(59)
  g <- 9L
  n_per <- 4L
  n <- g * n_per
  cluster <- rep(seq_len(g), each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  small_sample_warns <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_caught
  )
  msgs <- vapply(small_sample_warns, conditionMessage, character(1))
  # G=9 < 10: should say "highly unreliable" (strong warning)
  expect_true(any(grepl("highly unreliable", msgs, ignore.case = TRUE)))
  expect_true(any(grepl("G=9", msgs)))
})

# CV-06b: G=19 boundary triggers informational warning
test_that("CV-06b: G=19 triggers informational warning (10<=G<20)", {
  set.seed(619)
  g <- 19L
  n_per <- 4L
  n <- g * n_per
  cluster <- rep(seq_len(g), each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  small_sample_warns <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_caught
  )
  msgs <- vapply(small_sample_warns, conditionMessage, character(1))
  # G=19 in [10, 20): informational, not strong
  expect_false(any(grepl("highly unreliable", msgs, ignore.case = TRUE)))
  expect_true(any(grepl("Bootstrap", msgs, ignore.case = TRUE)))
  expect_true(any(grepl("G=19", msgs)))
})

# CV-08b: G=20 boundary (exactly 20) no cluster count warning
test_that("CV-08b: G=20 boundary emits no cluster count warning", {
  set.seed(820)
  g <- 20L
  n_per <- 3L
  n <- g * n_per
  cluster <- rep(seq_len(g), each = n_per)
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + rnorm(n, 0, 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, vce = "cluster", cluster = cluster),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # No "few clusters" or "Very few" warnings
  cluster_count_warns <- Filter(
    function(w) {
      msg <- conditionMessage(w)
      grepl("few cluster", msg, ignore.case = TRUE) ||
      grepl("Very few", msg, ignore.case = TRUE)
    },
    warnings_caught
  )
  expect_equal(length(cluster_count_warns), 0L)
})

# CV-03b: Single unit nesting violation detected
test_that("CV-03b: single unit nesting violation detected", {
  set.seed(32)
  n_units <- 10L
  t_periods <- 5L
  dt <- data.table::CJ(unit = 1:n_units, year = 1:t_periods)
  dt[, d := as.integer(unit <= 3)]
  dt[, post := as.integer(year >= 4)]
  dt[, y := rnorm(.N)]
  dt[, cluster := ifelse(unit <= 5, "A", "B")]
  # Only unit 3 violates nesting
  dt[unit == 3 & year >= 4, cluster := "B"]

  expect_error(
    suppressWarnings(
      lwdid(dt, y = "y", ivar = "unit", tvar = "year",
            d = "d", post = "post", rolling = "demean",
            vce = "cluster", cluster_var = "cluster")
    ),
    class = "lwdid_invalid_parameter"
  )
  err <- tryCatch(
    suppressWarnings(
      lwdid(dt, y = "y", ivar = "unit", tvar = "year",
            d = "d", post = "post", rolling = "demean",
            vce = "cluster", cluster_var = "cluster")
    ),
    error = function(e) e
  )
  # Should mention 1 unit violating
  expect_true(grepl("1 unit", conditionMessage(err)))
})

# CV-04b: Character cluster_var through lwdid() pipeline end-to-end
test_that("CV-04b: character cluster_var through lwdid() pipeline", {
  data(smoking)
  # state column is character/factor — use it directly
  result <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "state")
  )
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$vce_type, "cluster")
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se_att))
  expect_true(result$se_att > 0)
})
