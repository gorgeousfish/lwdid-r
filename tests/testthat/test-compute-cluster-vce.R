# ============================================================================
# Tests for compute_cluster_vce() — Story E3-03
# ============================================================================
# Test groups:
#   Group 1: sandwich::vcovCL consistency (E3-03.2)
#   (more groups will be appended by later tasks)
# ============================================================================

# Helper: suppress expected warnings from compute_cluster_vce
# (small sample warnings for G < 20)
suppress_cluster_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    lwdid_small_sample = function(w) invokeRestart("muffleWarning")
  )
}

# ============================================================================
# Group 1: sandwich::vcovCL consistency (E3-03.2)
# ============================================================================

#' Generate clustered test data for cluster-robust VCE tests
#'
#' @param n Total number of observations.
#' @param G Number of clusters (must divide n evenly).
#' @param seed Random seed.
#' @return data.frame with columns y, D, x, cl.
make_cluster_test_data <- function(n = 60, G = 10, seed = 42) {
  set.seed(seed)
  cluster <- rep(seq_len(G), each = n / G)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  # Cluster-level random effect for realistic within-cluster correlation
  cluster_effect <- rnorm(G, sd = 1.5)[cluster]
  y <- 1 + 2 * d + 0.5 * x + cluster_effect + rnorm(n)
  data.frame(y = y, D = d, x = x, cl = cluster)
}

test_that("cluster VCE matches sandwich::vcovCL(type='HC1') within 1e-12", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 42)
  fit <- lm(y ~ D + x, data = df)

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
  )
  expected <- sandwich::vcovCL(fit, cluster = df$cl, type = "HC1")

  expect_equal(
    result$vcov, expected,
    tolerance = 1e-12,
    label = "compute_cluster_vce vcov",
    expected.label = "sandwich::vcovCL(type='HC1')"
  )
})

test_that("df equals G - 1", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 42)
  fit <- lm(y ~ D + x, data = df)
  G <- length(unique(df$cl))

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
  )

  expect_identical(result$df, G - 1L)
})

test_that("n_clusters equals G", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 42)
  fit <- lm(y ~ D + x, data = df)
  G <- length(unique(df$cl))

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
  )

  expect_identical(result$n_clusters, G)
})

test_that("consistency across G values: G=5, 10, 20, 50", {
  for (g_val in c(5L, 10L, 20L, 50L)) {
    n <- g_val * 6L
    df <- make_cluster_test_data(n = n, G = g_val, seed = 42)
    fit <- lm(y ~ D + x, data = df)

    result <- suppress_cluster_warnings(
      compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
    )
    expected <- sandwich::vcovCL(
      fit, cluster = df$cl, type = "HC1"
    )

    # VCE matrix matches sandwich
    expect_equal(
      result$vcov, expected,
      tolerance = 1e-12,
      label = sprintf(
        "compute_cluster_vce vcov for G=%d", g_val
      ),
      expected.label = sprintf(
        "sandwich::vcovCL for G=%d", g_val
      )
    )

    # df = G - 1
    expect_identical(
      result$df, g_val - 1L,
      label = sprintf("df for G=%d", g_val)
    )

    # n_clusters = G
    expect_identical(
      result$n_clusters, g_val,
      label = sprintf("n_clusters for G=%d", g_val)
    )
  }
})

test_that("return structure: vcov is matrix, se derivable from vcov", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 42)
  fit <- lm(y ~ D + x, data = df)
  p <- length(coef(fit))

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
  )

  # Exactly three fields: vcov, df, n_clusters
  expect_identical(
    sort(names(result)),
    c("df", "n_clusters", "vcov")
  )

  # vcov is a p x p matrix
  expect_true(is.matrix(result$vcov))
  expect_equal(dim(result$vcov), c(p, p))

  # SE derivable from vcov: sqrt(diag(vcov)) should be positive
  se <- sqrt(diag(result$vcov))
  expect_true(is.numeric(se))
  expect_equal(length(se), p)
  expect_true(all(se > 0))

  # df is integer
  expect_true(is.integer(result$df) || is.numeric(result$df))

  # n_clusters is integer
  expect_true(
    is.integer(result$n_clusters) ||
      is.numeric(result$n_clusters)
  )
})

test_that("vcov matrix is symmetric and positive semi-definite", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 42)
  fit <- lm(y ~ D + x, data = df)

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
  )

  # Symmetric (sandwich::vcovCL may have tiny floating-point asymmetry)
  expect_equal(
    result$vcov, t(result$vcov),
    tolerance = 1e-12,
    label = "vcov symmetry"
  )

  # Positive semi-definite: all eigenvalues >= 0
  eig <- eigen(result$vcov, symmetric = TRUE)$values
  expect_true(
    all(eig >= -1e-10),
    info = paste(
      "min eigenvalue =", min(eig),
      "should be >= 0"
    )
  )
})

test_that("cluster VCE SE larger than OLS SE (with cluster correlation)", {
  # With cluster-level random effects, cluster-robust SE should
  # generally be larger than OLS SE
  df <- make_cluster_test_data(n = 60, G = 10, seed = 42)
  fit <- lm(y ~ D + x, data = df)

  se_ols <- sqrt(diag(vcov(fit)))
  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
  )
  se_cluster <- sqrt(diag(result$vcov))

  # At least the treatment effect SE should be larger
  # (cluster effects inflate within-cluster correlation)
  d_idx <- which(names(coef(fit)) == "D")
  expect_true(
    se_cluster[d_idx] > se_ols[d_idx],
    info = sprintf(
      "Cluster SE(D)=%.6f should > OLS SE(D)=%.6f",
      se_cluster[d_idx], se_ols[d_idx]
    )
  )
})

test_that("coefficient names preserved in vcov matrix", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 42)
  fit <- lm(y ~ D + x, data = df)

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
  )

  expect_identical(
    rownames(result$vcov),
    names(coef(fit))
  )
  expect_identical(
    colnames(result$vcov),
    names(coef(fit))
  )
})


# ============================================================================
# Group 2: Graded Cluster Warning Tests (E3-03.3)
# ============================================================================

test_that("G=5 (<10): triggers strong lwdid_small_sample warning", {
  df <- make_cluster_test_data(n = 30, G = 5, seed = 100)
  fit <- lm(y ~ D + x, data = df)

  expect_warning(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1"),
    regexp = "highly unreliable",
    class = "lwdid_small_sample"
  )
})

test_that("G=9 (<10 boundary): triggers strong lwdid_small_sample warning", {
  df <- make_cluster_test_data(n = 54, G = 9, seed = 101)
  fit <- lm(y ~ D + x, data = df)

  expect_warning(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1"),
    regexp = "highly unreliable",
    class = "lwdid_small_sample"
  )
})

test_that("G=10 (boundary): triggers informational warning, NOT strong", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 102)
  fit <- lm(y ~ D + x, data = df)

  # Should get informational warning (mentions "Wild Cluster Bootstrap")
  expect_warning(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1"),
    regexp = "Wild Cluster Bootstrap",
    class = "lwdid_small_sample"
  )

  # The message should NOT contain "highly unreliable"
  w <- tryCatch(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1"),
    lwdid_small_sample = function(w) w
  )
  expect_false(grepl("highly unreliable", w$message))
})

test_that("G=15 (10<=G<20): triggers informational lwdid_small_sample warning", {
  df <- make_cluster_test_data(n = 90, G = 15, seed = 103)
  fit <- lm(y ~ D + x, data = df)

  expect_warning(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1"),
    regexp = "Wild Cluster Bootstrap",
    class = "lwdid_small_sample"
  )
})

test_that("G=19 (<20 boundary): triggers informational lwdid_small_sample warning", {
  df <- make_cluster_test_data(n = 114, G = 19, seed = 104)
  fit <- lm(y ~ D + x, data = df)

  expect_warning(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1"),
    regexp = "Wild Cluster Bootstrap",
    class = "lwdid_small_sample"
  )
})

test_that("G=20 (>=20): NO lwdid_small_sample warning emitted", {
  df <- make_cluster_test_data(n = 120, G = 20, seed = 105)
  fit <- lm(y ~ D + x, data = df)

  # Use withCallingHandlers to fail if lwdid_small_sample is raised
  withCallingHandlers(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1"),
    lwdid_small_sample = function(w) {
      fail(sprintf(
        "Unexpected lwdid_small_sample warning for G=20: %s",
        conditionMessage(w)
      ))
    }
  )
  succeed()
})

# ============================================================================
# Group 3: Input Validation Tests (E3-03.4)
# ============================================================================

test_that("cluster with single NA throws lwdid_invalid_parameter", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 200)
  fit <- lm(y ~ D + x, data = df)
  cl_bad <- df$cl
  cl_bad[1] <- NA

  expect_error(
    compute_cluster_vce(fit, cluster = cl_bad, type = "HC1"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("cluster with multiple NAs throws lwdid_invalid_parameter", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 201)
  fit <- lm(y ~ D + x, data = df)
  cl_bad <- df$cl
  cl_bad[c(1, 5, 10)] <- NA

  expect_error(
    compute_cluster_vce(fit, cluster = cl_bad, type = "HC1"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("cluster all NA throws lwdid_invalid_parameter", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 202)
  fit <- lm(y ~ D + x, data = df)
  cl_bad <- rep(NA, nrow(df))

  expect_error(
    compute_cluster_vce(fit, cluster = cl_bad, type = "HC1"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("cluster length too long throws lwdid_invalid_parameter", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 203)
  fit <- lm(y ~ D + x, data = df)
  cl_long <- c(df$cl, 1L, 2L, 3L)

  expect_error(
    compute_cluster_vce(fit, cluster = cl_long, type = "HC1"),
    regexp = "length",
    class = "lwdid_invalid_parameter"
  )
})

test_that("cluster length too short throws lwdid_invalid_parameter", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 204)
  fit <- lm(y ~ D + x, data = df)
  cl_short <- df$cl[1:30]

  expect_error(
    compute_cluster_vce(fit, cluster = cl_short, type = "HC1"),
    regexp = "length",
    class = "lwdid_invalid_parameter"
  )
})

test_that("cluster empty vector throws lwdid_invalid_parameter (length mismatch)", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 205)
  fit <- lm(y ~ D + x, data = df)

  expect_error(
    compute_cluster_vce(fit, cluster = integer(0), type = "HC1"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("G=1 (single cluster) throws lwdid_insufficient_data", {
  set.seed(300)
  n <- 20
  y <- rnorm(n)
  x <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  fit <- lm(y ~ d + x)
  cl_single <- rep(1L, n)

  expect_error(
    compute_cluster_vce(fit, cluster = cl_single, type = "HC1"),
    class = "lwdid_insufficient_data"
  )
})

test_that("factor with unused levels: G based on actual values, G=1 throws error", {
  set.seed(301)
  n <- 20
  y <- rnorm(n)
  x <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  fit <- lm(y ~ d + x)
  # Factor with 3 levels but only "A" is used -> unique() gives 1 actual cluster
  cl_factor <- factor(rep("A", n), levels = c("A", "B", "C"))

  expect_error(
    compute_cluster_vce(fit, cluster = cl_factor, type = "HC1"),
    class = "lwdid_insufficient_data"
  )
})

test_that("G=2 boundary: normal computation, df=1, n_clusters=2", {
  set.seed(302)
  n <- 20
  y <- rnorm(n)
  x <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  fit <- lm(y ~ d + x)
  cl_two <- rep(1:2, each = n / 2)

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = cl_two, type = "HC1")
  )

  expect_equal(result$df, 1L)
  expect_equal(result$n_clusters, 2L)
  expect_true(is.matrix(result$vcov))
})

test_that("validation order: NA error comes before length mismatch", {
  df <- make_cluster_test_data(n = 60, G = 10, seed = 206)
  fit <- lm(y ~ D + x, data = df)
  # Cluster with both NA and wrong length (too long, with NAs)
  cl_bad <- c(df$cl, NA, NA, NA)

  # Should throw lwdid_invalid_parameter for NA (checked first)
  expect_error(
    compute_cluster_vce(fit, cluster = cl_bad, type = "HC1"),
    regexp = "NA",
    class = "lwdid_invalid_parameter"
  )
})


# ============================================================================
# Group 4: Cluster Size Imbalance Tests (E3-03.5)
# ============================================================================

test_that("high imbalance (CV >> 1): triggers imbalance warning", {
  # G=25 (>=20, no graded small-sample warning)

  # 1 cluster with 49 obs + 24 clusters with 1 obs each => n=73
  set.seed(400)
  n <- 73L
  cl <- c(rep(1L, 49L), 2:25)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  fit <- lm(y ~ d + x)

  # Verify CV is indeed > 1.0
  sizes <- as.integer(table(cl))
  cv <- sd(sizes) / mean(sizes)
  expect_true(cv > 1.0, info = sprintf("CV=%.4f should be > 1.0", cv))

  # Should trigger lwdid_small_sample warning with "unbalanced" or
  # "imbalance" in the message
  warnings_caught <- list()
  withCallingHandlers(
    compute_cluster_vce(fit, cluster = cl, type = "HC1"),
    lwdid_small_sample = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  # At least one warning should mention imbalance
  imbalance_msgs <- vapply(warnings_caught, function(w) {
    grepl("unbalanced|imbalance", conditionMessage(w),
          ignore.case = TRUE)
  }, logical(1))
  expect_true(
    any(imbalance_msgs),
    info = "Expected at least one warning about cluster imbalance"
  )
})

test_that("balanced clusters (CV=0): no imbalance warning", {
  # G=20 (>=20, no graded warning), each cluster 3 obs => n=60
  df <- make_cluster_test_data(n = 60, G = 20, seed = 401)
  fit <- lm(y ~ D + x, data = df)

  # Verify CV is 0
  sizes <- as.integer(table(df$cl))
  cv <- sd(sizes) / mean(sizes)
  expect_equal(cv, 0, info = "Balanced clusters should have CV=0")

  # No warnings at all for G>=20 with balanced clusters
  withCallingHandlers(
    compute_cluster_vce(fit, cluster = df$cl, type = "HC1"),
    lwdid_small_sample = function(w) {
      fail(sprintf(
        "Unexpected lwdid_small_sample warning: %s",
        conditionMessage(w)
      ))
    }
  )
  succeed()
})

test_that("boundary: CV just above 1.0 triggers, just below does not", {
  # --- CV > 1.0: should trigger ---
  # G=20, construct sizes so CV slightly > 1.0
  # 1 cluster with 40 obs + 19 clusters with 2 obs => n=78
  set.seed(402)
  n_above <- 78L
  cl_above <- c(rep(1L, 40L), rep(2:20, each = 2L))
  sizes_above <- as.integer(table(cl_above))
  cv_above <- sd(sizes_above) / mean(sizes_above)
  expect_true(cv_above > 1.0,
    info = sprintf("CV_above=%.4f should be > 1.0", cv_above))

  d <- rbinom(n_above, 1, 0.5)
  x <- rnorm(n_above)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n_above)
  fit_above <- lm(y ~ d + x)

  got_imbalance_above <- FALSE
  withCallingHandlers(
    compute_cluster_vce(fit_above, cluster = cl_above, type = "HC1"),
    lwdid_small_sample = function(w) {
      if (grepl("unbalanced|imbalance", conditionMessage(w),
                ignore.case = TRUE)) {
        got_imbalance_above <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_true(got_imbalance_above,
    info = "CV > 1.0 should trigger imbalance warning")


  # --- CV <= 1.0: should NOT trigger ---
  # G=20, each cluster 3 obs => n=60, CV=0
  set.seed(403)
  n_below <- 60L
  cl_below <- rep(1:20, each = 3L)
  sizes_below <- as.integer(table(cl_below))
  cv_below <- sd(sizes_below) / mean(sizes_below)
  expect_true(cv_below <= 1.0,
    info = sprintf("CV_below=%.4f should be <= 1.0", cv_below))

  d2 <- rbinom(n_below, 1, 0.5)
  x2 <- rnorm(n_below)
  y2 <- 1 + 2 * d2 + 0.5 * x2 + rnorm(n_below)
  fit_below <- lm(y2 ~ d2 + x2)

  withCallingHandlers(
    compute_cluster_vce(fit_below, cluster = cl_below, type = "HC1"),
    lwdid_small_sample = function(w) {
      if (grepl("unbalanced|imbalance", conditionMessage(w),
                ignore.case = TRUE)) {
        fail("CV <= 1.0 should NOT trigger imbalance warning")
      }
      invokeRestart("muffleWarning")
    }
  )
  succeed()
})

test_that("single-obs clusters (G=N, CV=0): no imbalance, computes OK", {
  # Each observation is its own cluster => G=N=30, CV=0
  set.seed(404)
  n <- 30L
  cl <- seq_len(n)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  fit <- lm(y ~ d + x)

  # Verify CV = 0
  sizes <- as.integer(table(cl))
  cv <- sd(sizes) / mean(sizes)
  expect_equal(cv, 0, info = "Single-obs clusters should have CV=0")

  # No imbalance warning (G=30 >= 20, so no graded warning either)
  withCallingHandlers(
    compute_cluster_vce(fit, cluster = cl, type = "HC1"),
    lwdid_small_sample = function(w) {
      fail(sprintf(
        "Unexpected lwdid_small_sample warning: %s",
        conditionMessage(w)
      ))
    }
  )

  # Verify result structure and correctness
  result <- compute_cluster_vce(fit, cluster = cl, type = "HC1")
  expect_identical(result$n_clusters, length(unique(cl)))
  expect_identical(result$df, n - 1L)
  expect_true(is.matrix(result$vcov))
  expect_equal(dim(result$vcov), c(3L, 3L))

  # Cross-check with sandwich
  expected <- sandwich::vcovCL(fit, cluster = cl, type = "HC1")
  expect_equal(result$vcov, expected, tolerance = 1e-12)
})


# ============================================================================
# Group 5: Character/Numeric/Factor Cluster Variable Tests (E3-03.6)
# ============================================================================

test_that("character cluster variable: correct computation", {
  set.seed(500)
  states <- c("CA", "NY", "TX", "FL", "WA")
  n_per <- 12L
  n <- length(states) * n_per
  cl <- rep(states, each = n_per)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  cluster_effect <- rnorm(length(states), sd = 1.5)
  names(cluster_effect) <- states
  y <- 1 + 2 * d + 0.5 * x + cluster_effect[cl] + rnorm(n)
  fit <- lm(y ~ d + x)

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = cl, type = "HC1")
  )

  expect_identical(result$n_clusters, length(states))
  expect_identical(result$df, length(states) - 1L)
  expect_true(is.matrix(result$vcov))
  expect_equal(dim(result$vcov), c(3L, 3L))

  # Cross-check with sandwich
  expected <- sandwich::vcovCL(fit, cluster = cl, type = "HC1")
  expect_equal(result$vcov, expected, tolerance = 1e-12)
})

test_that("numeric cluster variable: correct computation", {
  set.seed(501)
  g_ids <- 1:5
  n_per <- 12L
  n <- length(g_ids) * n_per
  cl <- rep(g_ids, each = n_per)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  cluster_effect <- rnorm(length(g_ids), sd = 1.5)
  y <- 1 + 2 * d + 0.5 * x + cluster_effect[cl] + rnorm(n)
  fit <- lm(y ~ d + x)

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = cl, type = "HC1")
  )

  expect_identical(result$n_clusters, length(g_ids))
  expect_identical(result$df, length(g_ids) - 1L)
  expect_true(is.matrix(result$vcov))

  expected <- sandwich::vcovCL(fit, cluster = cl, type = "HC1")
  expect_equal(result$vcov, expected, tolerance = 1e-12)
})

test_that("factor cluster variable: correct computation", {
  set.seed(502)
  g_labels <- c("alpha", "beta", "gamma", "delta", "epsilon")
  n_per <- 12L
  n <- length(g_labels) * n_per
  cl <- factor(rep(g_labels, each = n_per), levels = g_labels)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  cluster_effect <- rnorm(length(g_labels), sd = 1.5)
  names(cluster_effect) <- g_labels
  y <- 1 + 2 * d + 0.5 * x + cluster_effect[as.character(cl)] +
    rnorm(n)
  fit <- lm(y ~ d + x)

  result <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = cl, type = "HC1")
  )

  expect_identical(result$n_clusters, length(g_labels))
  expect_identical(result$df, length(g_labels) - 1L)
  expect_true(is.matrix(result$vcov))

  expected <- sandwich::vcovCL(fit, cluster = cl, type = "HC1")
  expect_equal(result$vcov, expected, tolerance = 1e-12)
})

test_that("character vs numeric cluster: identical VCE matrix", {
  set.seed(503)
  n_per <- 12L
  g <- 5L
  n <- g * n_per

  # Same data for both
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  cluster_effect <- rnorm(g, sd = 1.5)
  y <- 1 + 2 * d + 0.5 * x +
    cluster_effect[rep(seq_len(g), each = n_per)] + rnorm(n)
  fit <- lm(y ~ d + x)

  # Numeric cluster: 1, 2, 3, 4, 5
  cl_num <- rep(seq_len(g), each = n_per)
  result_num <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = cl_num, type = "HC1")
  )

  # Character cluster: "1", "2", "3", "4", "5"
  cl_chr <- as.character(cl_num)
  result_chr <- suppress_cluster_warnings(
    compute_cluster_vce(fit, cluster = cl_chr, type = "HC1")
  )

  # VCE matrices should be identical (< 1e-12)
  expect_equal(
    result_num$vcov, result_chr$vcov,
    tolerance = 1e-12,
    label = "numeric cluster VCE",
    expected.label = "character cluster VCE"
  )

  # Metadata should also match
  expect_identical(result_num$df, result_chr$df)
  expect_identical(result_num$n_clusters, result_chr$n_clusters)
})


# ============================================================================
# Group 6: Numerical Verification Tests (E3-03.7)
# ============================================================================

# --- 6.1 Cluster SE > Homoskedastic SE with strong within-cluster correlation ---

test_that("cluster SE > OLS SE for all coefficients with cluster-level regressors", {
  # G=20 balanced clusters with strong within-cluster correlation.
  # Treatment and covariate assigned at CLUSTER level so that
  # within-cluster correlation inflates SE for ALL coefficients.
  set.seed(600)
  G <- 20L
  n_per <- 15L
  n <- G * n_per
  cl <- rep(seq_len(G), each = n_per)

  # Cluster-level treatment and covariate (constant within cluster)
  d_cl <- rbinom(G, 1, 0.5)
  x_cl <- rnorm(G)
  d <- d_cl[cl]
  x <- x_cl[cl]

  # Strong cluster effect (sd=4) creates high within-cluster correlation
  cluster_effect <- rnorm(G, sd = 4)[cl]
  y <- 1 + 2 * d + 0.5 * x + cluster_effect + rnorm(n, sd = 0.5)
  fit <- lm(y ~ d + x)

  se_ols <- sqrt(diag(vcov(fit)))
  result <- compute_cluster_vce(fit, cluster = cl, type = "HC1")
  se_cluster <- sqrt(diag(result$vcov))

  # ALL coefficients should have cluster SE > OLS SE
  for (i in seq_along(se_ols)) {
    expect_true(
      se_cluster[i] > se_ols[i],
      info = sprintf(
        "Coef '%s': cluster SE=%.6f should > OLS SE=%.6f",
        names(coef(fit))[i], se_cluster[i], se_ols[i]
      )
    )
  }
})

test_that("cluster SE ratio > 1 increases with cluster effect strength", {
  # Verify that stronger within-cluster correlation leads to larger
  # cluster-to-OLS SE ratio for the intercept (cluster-level quantity).
  set.seed(601)
  G <- 20L
  n_per <- 15L
  n <- G * n_per
  cl <- rep(seq_len(G), each = n_per)

  # Cluster-level regressors for reliable SE inflation
  d_cl <- rbinom(G, 1, 0.5)
  x_cl <- rnorm(G)
  d <- d_cl[cl]
  x <- x_cl[cl]

  ratios <- numeric(3)
  # Widely separated sd values for robust monotonicity
  sd_vals <- c(0.5, 3, 10)

  for (k in seq_along(sd_vals)) {
    set.seed(601)
    cluster_effect <- rnorm(G, sd = sd_vals[k])[cl]
    y <- 1 + 2 * d + 0.5 * x + cluster_effect + rnorm(n, sd = 1)
    fit <- lm(y ~ d + x)

    se_ols <- sqrt(diag(vcov(fit)))
    result <- compute_cluster_vce(fit, cluster = cl, type = "HC1")
    se_cluster <- sqrt(diag(result$vcov))

    # Use intercept ratio as representative
    ratios[k] <- se_cluster[1] / se_ols[1]
  }

  # All ratios should be > 1
  for (k in seq_along(ratios)) {
    expect_true(ratios[k] > 1,
      info = sprintf("ratio(sd=%.1f)=%.4f should > 1",
                     sd_vals[k], ratios[k]))
  }

  # Ratios should be monotonically increasing
  expect_true(ratios[2] > ratios[1],
    info = sprintf("ratio(sd=3)=%.4f should > ratio(sd=0.5)=%.4f",
                   ratios[2], ratios[1]))
  expect_true(ratios[3] > ratios[2],
    info = sprintf("ratio(sd=10)=%.4f should > ratio(sd=3)=%.4f",
                   ratios[3], ratios[2]))
})

# --- 6.2 Small sample correction factor verification ---

test_that("HC1 vcov = HC0(cadjust=FALSE) * G/(G-1) * (N-1)/(N-p)", {
  # Use G=20 balanced clusters to avoid small sample warnings.
  # Verify the correction factor relationship between HC1 and HC0.
  set.seed(610)
  G <- 20L
  n_per <- 10L
  n <- G * n_per
  cl <- rep(seq_len(G), each = n_per)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  cluster_effect <- rnorm(G, sd = 2)[cl]
  y <- 1 + 2 * d + 0.5 * x + cluster_effect + rnorm(n)
  fit <- lm(y ~ d + x)
  p <- length(coef(fit))

  # HC1 from compute_cluster_vce (our function)
  result_hc1 <- compute_cluster_vce(fit, cluster = cl, type = "HC1")

  # HC0 from sandwich::vcovCL with cadjust=FALSE
  # (cadjust=FALSE disables the G/(G-1) correction in sandwich)
  vcov_hc0_raw <- sandwich::vcovCL(
    fit, cluster = cl, type = "HC0", cadjust = FALSE
  )

  # The correction factor for HC1 with cadjust=TRUE (default):
  # G/(G-1) * (N-1)/(N-p)
  correction <- (G / (G - 1)) * ((n - 1) / (n - p))

  # vcov_HC1 should equal vcov_HC0_raw * correction
  expected_hc1 <- vcov_hc0_raw * correction
  expect_equal(
    result_hc1$vcov, expected_hc1,
    tolerance = 1e-10,
    label = "HC1 vcov via correction factor",
    expected.label = "HC0(cadjust=FALSE) * G/(G-1) * (N-1)/(N-p)"
  )
})

test_that("correction factor numerical value for known G=20, N=200, p=3", {
  # Verify the correction factor itself for a specific case
  G <- 20L
  n <- 200L
  p <- 3L  # intercept + d + x

  correction <- (G / (G - 1)) * ((n - 1) / (n - p))

  # Manual calculation: 20/19 * 199/197
  expected_correction <- (20 / 19) * (199 / 197)
  expect_equal(correction, expected_correction, tolerance = 1e-14)

  # Verify it's > 1 (correction inflates variance)
  expect_true(correction > 1,
    info = sprintf("Correction=%.6f should be > 1", correction))

  # Verify numerical value with vibe-math cross-check:
  # 20/19 = 1.052631..., 199/197 = 1.010152..., product ~ 1.063...
  expect_true(correction > 1.06 && correction < 1.07,
    info = sprintf("Correction=%.6f should be ~1.063", correction))
})

test_that("compute_cluster_vce HC1 matches sandwich::vcovCL exactly", {
  # Direct comparison for multiple configurations
  configs <- list(
    list(G = 20L, n_per = 5L,  seed = 620),
    list(G = 30L, n_per = 10L, seed = 621),
    list(G = 50L, n_per = 8L,  seed = 622)
  )

  for (cfg in configs) {
    n <- cfg$G * cfg$n_per
    df <- make_cluster_test_data(n = n, G = cfg$G, seed = cfg$seed)
    fit <- lm(y ~ D + x, data = df)

    result <- suppress_cluster_warnings(
      compute_cluster_vce(fit, cluster = df$cl, type = "HC1")
    )
    expected <- sandwich::vcovCL(
      fit, cluster = df$cl, type = "HC1"
    )

    expect_equal(
      result$vcov, expected,
      tolerance = 1e-12,
      label = sprintf("vcov for G=%d, n=%d", cfg$G, n),
      expected.label = "sandwich::vcovCL"
    )
  }
})

# --- 6.3 Smoking dataset cluster SE ---

test_that("smoking data: cluster SE matches sandwich::vcovCL (< 1e-12)", {
  data(smoking, package = "lwdid")

  # Regression: cigsale ~ d + retprice, clustered by state
  sm <- smoking[complete.cases(
    smoking[, c("cigsale", "d", "retprice", "state")]
  ), ]
  fit <- lm(cigsale ~ d + retprice, data = sm)

  result <- compute_cluster_vce(
    fit, cluster = sm$state, type = "HC1"
  )
  expected <- sandwich::vcovCL(
    fit, cluster = sm$state, type = "HC1"
  )

  expect_equal(
    result$vcov, expected,
    tolerance = 1e-12,
    label = "smoking cluster vcov",
    expected.label = "sandwich::vcovCL on smoking data"
  )
})

test_that("smoking data: df = G - 1 where G = number of unique states", {
  data(smoking, package = "lwdid")

  sm <- smoking[complete.cases(
    smoking[, c("cigsale", "d", "retprice", "state")]
  ), ]
  fit <- lm(cigsale ~ d + retprice, data = sm)
  G <- length(unique(sm$state))

  result <- compute_cluster_vce(
    fit, cluster = sm$state, type = "HC1"
  )

  expect_identical(result$df, G - 1L)
  expect_identical(result$n_clusters, G)
})

test_that("smoking data: cluster SE > OLS SE for intercept and retprice", {
  # Note: the treatment variable 'd' varies within states (time-varying),

  # so cluster SE for 'd' can be smaller than OLS SE. However, the
  # intercept and retprice (which have strong between-state variation)
  # should have cluster SE > OLS SE.
  data(smoking, package = "lwdid")

  sm <- smoking[complete.cases(
    smoking[, c("cigsale", "d", "retprice", "state")]
  ), ]
  fit <- lm(cigsale ~ d + retprice, data = sm)

  se_ols <- sqrt(diag(vcov(fit)))
  result <- compute_cluster_vce(
    fit, cluster = sm$state, type = "HC1"
  )
  se_cluster <- sqrt(diag(result$vcov))

  # Intercept: strong between-state variation
  intercept_idx <- which(names(coef(fit)) == "(Intercept)")
  expect_true(
    se_cluster[intercept_idx] > se_ols[intercept_idx],
    info = sprintf(
      "Intercept: cluster SE=%.4f should > OLS SE=%.4f",
      se_cluster[intercept_idx], se_ols[intercept_idx]
    )
  )

  # retprice: between-state price variation
  rp_idx <- which(names(coef(fit)) == "retprice")
  expect_true(
    se_cluster[rp_idx] > se_ols[rp_idx],
    info = sprintf(
      "retprice: cluster SE=%.4f should > OLS SE=%.4f",
      se_cluster[rp_idx], se_ols[rp_idx]
    )
  )
})

test_that("smoking data: cluster SEs are numerically reasonable", {
  data(smoking, package = "lwdid")

  sm <- smoking[complete.cases(
    smoking[, c("cigsale", "d", "retprice", "state")]
  ), ]
  fit <- lm(cigsale ~ d + retprice, data = sm)

  result <- compute_cluster_vce(
    fit, cluster = sm$state, type = "HC1"
  )
  se_cluster <- sqrt(diag(result$vcov))

  # All SEs should be positive and finite
  expect_true(all(se_cluster > 0))
  expect_true(all(is.finite(se_cluster)))

  # SEs should be reasonable (cigsale ~ 100 packs per capita)
  expect_true(all(se_cluster < 1000),
    info = "Cluster SEs should be numerically reasonable")

  # Placeholder: Stata benchmark values for future cross-validation
  # Stata: regress cigsale d retprice, vce(cluster state)
  # TODO: Add Stata benchmark SE values when available
  # stata_se_intercept <- ...
  # stata_se_d <- ...
  # stata_se_retprice <- ...
})


