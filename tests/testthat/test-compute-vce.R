# ============================================================================
# Tests for compute_vce()
# Story E3-01: Variance-Covariance Estimation Dispatcher
#
# Test Groups:
#   Group 1: VCE type dispatch correctness (7 HC + 1 cluster)
#   Group 2: Return value structure
#   Group 3: Case-insensitivity
# ============================================================================

# ============================================================================
# Helper functions
# ============================================================================

#' Generate simple test data for VCE dispatch tests
make_test_data <- function(n = 50, seed = 42) {
  set.seed(seed)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  data.frame(y = y, D = d, x = x)
}

#' Generate clustered test data for cluster-robust VCE tests
make_cluster_data <- function(n = 100, G = 20, seed = 42) {
  set.seed(seed)
  cluster <- rep(1:G, each = n / G)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  data.frame(y = y, D = d, x = x, cl = cluster)
}

# ============================================================================
# Group 1: VCE type dispatch correctness
# ============================================================================

test_that("vce=NULL dispatches to homoskedastic OLS vcov", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = NULL)
  expected_vcov <- vcov(fit)
  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)
  expect_equal(res$df, fit$df.residual)
  expect_identical(res$vce_type, "homoskedastic")
})

test_that("vce='hc0' dispatches to sandwich HC0", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "hc0")
  expected_vcov <- sandwich::vcovHC(fit, type = "HC0")
  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)
  expect_equal(res$df, fit$df.residual)
  expect_identical(res$vce_type, "HC0")
})

test_that("vce='hc1' dispatches to sandwich HC1", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "hc1")
  expected_vcov <- sandwich::vcovHC(fit, type = "HC1")
  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)
  expect_equal(res$df, fit$df.residual)
  expect_identical(res$vce_type, "HC1")
})

test_that("vce='robust' is alias for HC1", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "robust")
  expected_vcov <- sandwich::vcovHC(fit, type = "HC1")
  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)
  expect_equal(res$df, fit$df.residual)
  expect_identical(res$vce_type, "HC1")
})

test_that("vce='hc2' dispatches to sandwich HC2", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "hc2")
  expected_vcov <- sandwich::vcovHC(fit, type = "HC2")
  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)
  expect_equal(res$df, fit$df.residual)
  expect_identical(res$vce_type, "HC2")
})

test_that("vce='hc3' dispatches to sandwich HC3", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "hc3")
  expected_vcov <- sandwich::vcovHC(fit, type = "HC3")
  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)
  expect_equal(res$df, fit$df.residual)
  expect_identical(res$vce_type, "HC3")
})

test_that("vce='hc4' dispatches to sandwich HC4", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "hc4")
  expected_vcov <- sandwich::vcovHC(fit, type = "HC4")
  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)
  expect_equal(res$df, fit$df.residual)
  expect_identical(res$vce_type, "HC4")
})

test_that("vce='cluster' dispatches to sandwich vcovCL with HC1", {
  df <- make_cluster_data(n = 100, G = 20)
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "cluster", cluster = df$cl)
  expected_vcov <- sandwich::vcovCL(fit, cluster = df$cl, type = "HC1")
  G <- length(unique(df$cl))
  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)
  expect_equal(res$df, G - 1L)
  expect_identical(res$vce_type, "cluster")
  expect_identical(res$n_clusters, G)
})

# ============================================================================
# Group 2: Return value structure
# ============================================================================

test_that("vcov matrix dimensions equal p x p", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  p <- length(coef(fit))

  # Homoskedastic
  res_homo <- compute_vce(fit, vce = NULL)
  expect_equal(dim(res_homo$vcov), c(p, p))

  # HC3
  res_hc3 <- compute_vce(fit, vce = "hc3")
  expect_equal(dim(res_hc3$vcov), c(p, p))

  # Cluster
  df_cl <- make_cluster_data(n = 100, G = 20)
  fit_cl <- lm(y ~ D + x, data = df_cl)
  p_cl <- length(coef(fit_cl))
  res_cl <- compute_vce(fit_cl, vce = "cluster", cluster = df_cl$cl)
  expect_equal(dim(res_cl$vcov), c(p_cl, p_cl))
})

test_that("vce_type field is correct string for each path", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)

  expect_identical(compute_vce(fit, vce = NULL)$vce_type,
                   "homoskedastic")
  expect_identical(compute_vce(fit, vce = "hc0")$vce_type, "HC0")
  expect_identical(compute_vce(fit, vce = "hc1")$vce_type, "HC1")
  expect_identical(compute_vce(fit, vce = "robust")$vce_type, "HC1")
  expect_identical(compute_vce(fit, vce = "hc2")$vce_type, "HC2")
  expect_identical(compute_vce(fit, vce = "hc3")$vce_type, "HC3")
  expect_identical(compute_vce(fit, vce = "hc4")$vce_type, "HC4")

  df_cl <- make_cluster_data(n = 100, G = 20)
  fit_cl <- lm(y ~ D + x, data = df_cl)
  expect_identical(
    compute_vce(fit_cl, vce = "cluster", cluster = df_cl$cl)$vce_type,
    "cluster"
  )
})

test_that("n_clusters is G for cluster VCE, NULL for others", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)

  # Non-cluster paths return NULL
  expect_null(compute_vce(fit, vce = NULL)$n_clusters)
  expect_null(compute_vce(fit, vce = "hc0")$n_clusters)
  expect_null(compute_vce(fit, vce = "hc1")$n_clusters)
  expect_null(compute_vce(fit, vce = "robust")$n_clusters)
  expect_null(compute_vce(fit, vce = "hc2")$n_clusters)
  expect_null(compute_vce(fit, vce = "hc3")$n_clusters)
  expect_null(compute_vce(fit, vce = "hc4")$n_clusters)

  # Cluster path returns G
  df_cl <- make_cluster_data(n = 100, G = 20)
  fit_cl <- lm(y ~ D + x, data = df_cl)
  G <- length(unique(df_cl$cl))
  res_cl <- compute_vce(fit_cl, vce = "cluster", cluster = df_cl$cl)
  expect_identical(res_cl$n_clusters, G)
})

test_that("df is fit$df.residual for non-cluster, G-1 for cluster", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)

  # Non-cluster paths use fit$df.residual
  expect_equal(compute_vce(fit, vce = NULL)$df, fit$df.residual)
  expect_equal(compute_vce(fit, vce = "hc0")$df, fit$df.residual)
  expect_equal(compute_vce(fit, vce = "hc1")$df, fit$df.residual)
  expect_equal(compute_vce(fit, vce = "robust")$df, fit$df.residual)
  expect_equal(compute_vce(fit, vce = "hc2")$df, fit$df.residual)
  expect_equal(compute_vce(fit, vce = "hc3")$df, fit$df.residual)
  expect_equal(compute_vce(fit, vce = "hc4")$df, fit$df.residual)

  # Cluster path uses G - 1
  df_cl <- make_cluster_data(n = 100, G = 20)
  fit_cl <- lm(y ~ D + x, data = df_cl)
  G <- length(unique(df_cl$cl))
  res_cl <- compute_vce(fit_cl, vce = "cluster", cluster = df_cl$cl)
  expect_equal(res_cl$df, G - 1L)
})

# ============================================================================
# Group 3: Case-insensitivity
# ============================================================================

test_that("HC3/hc3/Hc3 all dispatch to HC3 (same vcov)", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res_lower <- compute_vce(fit, vce = "hc3")
  res_upper <- compute_vce(fit, vce = "HC3")
  res_mixed <- compute_vce(fit, vce = "Hc3")
  expect_equal(res_lower$vcov, res_upper$vcov, tolerance = 1e-12)
  expect_equal(res_lower$vcov, res_mixed$vcov, tolerance = 1e-12)
  expect_identical(res_lower$vce_type, "HC3")
  expect_identical(res_upper$vce_type, "HC3")
  expect_identical(res_mixed$vce_type, "HC3")
})

test_that("ROBUST/Robust/robust all dispatch to HC1 (same vcov)", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  res_lower <- compute_vce(fit, vce = "robust")
  res_upper <- compute_vce(fit, vce = "ROBUST")
  res_mixed <- compute_vce(fit, vce = "Robust")
  expect_equal(res_lower$vcov, res_upper$vcov, tolerance = 1e-12)
  expect_equal(res_lower$vcov, res_mixed$vcov, tolerance = 1e-12)
  expect_identical(res_lower$vce_type, "HC1")
  expect_identical(res_upper$vce_type, "HC1")
  expect_identical(res_mixed$vce_type, "HC1")
})

test_that("CLUSTER/Cluster/cluster all dispatch to cluster VCE", {
  df_cl <- make_cluster_data(n = 100, G = 20)
  fit_cl <- lm(y ~ D + x, data = df_cl)
  res_lower <- compute_vce(fit_cl, vce = "cluster", cluster = df_cl$cl)
  res_upper <- compute_vce(fit_cl, vce = "CLUSTER", cluster = df_cl$cl)
  res_mixed <- compute_vce(fit_cl, vce = "Cluster", cluster = df_cl$cl)
  expect_equal(res_lower$vcov, res_upper$vcov, tolerance = 1e-12)
  expect_equal(res_lower$vcov, res_mixed$vcov, tolerance = 1e-12)
  expect_identical(res_lower$vce_type, "cluster")
  expect_identical(res_upper$vce_type, "cluster")
  expect_identical(res_mixed$vce_type, "cluster")
})

# ============================================================================
# Group 4: Error handling
# ============================================================================

test_that("invalid VCE type 'hc5' raises lwdid_invalid_vce error", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  expect_error(
    compute_vce(fit, vce = "hc5"),
    class = "lwdid_invalid_vce"
  )
})

test_that("invalid VCE type 'invalid' raises lwdid_invalid_vce error", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  expect_error(
    compute_vce(fit, vce = "invalid"),
    class = "lwdid_invalid_vce"
  )
})

test_that("empty string VCE type raises lwdid_invalid_vce error", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  expect_error(
    compute_vce(fit, vce = ""),
    class = "lwdid_invalid_vce"
  )
})

test_that("vce='cluster' with cluster=NULL raises lwdid_invalid_parameter error", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  expect_error(
    compute_vce(fit, vce = "cluster", cluster = NULL),
    class = "lwdid_invalid_parameter"
  )
})

test_that("cluster containing NA raises lwdid_invalid_parameter error", {
  df <- make_cluster_data(n = 100, G = 20)
  fit <- lm(y ~ D + x, data = df)
  cl_with_na <- df$cl
  cl_with_na[1] <- NA
  expect_error(
    compute_vce(fit, vce = "cluster", cluster = cl_with_na),
    class = "lwdid_invalid_parameter"
  )
})

test_that("cluster length mismatch raises lwdid_invalid_parameter error", {
  df <- make_cluster_data(n = 100, G = 20)
  fit <- lm(y ~ D + x, data = df)
  expect_error(
    compute_vce(fit, vce = "cluster", cluster = df$cl[1:50]),
    class = "lwdid_invalid_parameter"
  )
})

test_that("G < 2 (single cluster) raises lwdid_insufficient_data error", {
  df <- make_test_data(n = 30)
  fit <- lm(y ~ D + x, data = df)
  single_cluster <- rep(1L, 30)
  expect_error(
    compute_vce(fit, vce = "cluster", cluster = single_cluster),
    class = "lwdid_insufficient_data"
  )
})

# ============================================================================
# Group 5: Cluster graded warnings
# ============================================================================

#' Generate clustered data with a specific number of clusters G
make_cluster_data_g <- function(G, n_per_cluster = 10, seed = 42) {
  set.seed(seed)
  n <- G * n_per_cluster
  cluster <- rep(1:G, each = n_per_cluster)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  data.frame(y = y, D = d, x = x, cl = cluster)
}

test_that("G < 10 (G=5) triggers strong lwdid_small_sample warning with 'highly unreliable'", {
  df <- make_cluster_data_g(G = 5)
  fit <- lm(y ~ D + x, data = df)
  expect_warning(
    compute_vce(fit, vce = "cluster", cluster = df$cl),
    regexp = "highly unreliable",
    class = "lwdid_small_sample"
  )
})

test_that("10 <= G < 20 (G=15) triggers informational lwdid_small_sample warning with 'Wild Cluster Bootstrap'", {
  df <- make_cluster_data_g(G = 15)
  fit <- lm(y ~ D + x, data = df)
  expect_warning(
    compute_vce(fit, vce = "cluster", cluster = df$cl),
    regexp = "Wild Cluster Bootstrap",
    class = "lwdid_small_sample"
  )
})

test_that("G >= 20 does not trigger lwdid_small_sample warning", {
  df <- make_cluster_data_g(G = 25)
  fit <- lm(y ~ D + x, data = df)
  expect_no_warning(
    compute_vce(fit, vce = "cluster", cluster = df$cl)
  )
})

test_that("G = 10 boundary triggers informational warning (not strong), matches 'Wild Cluster Bootstrap'", {
  df <- make_cluster_data_g(G = 10)
  fit <- lm(y ~ D + x, data = df)
  # G=10 is in [10, 20) range -> informational, not strong
  w <- tryCatch(
    compute_vce(fit, vce = "cluster", cluster = df$cl),
    lwdid_small_sample = function(w) w
  )
  expect_s3_class(w, "lwdid_small_sample")
  expect_match(conditionMessage(w), "Wild Cluster Bootstrap")
  # Should NOT contain "highly unreliable"
  expect_false(grepl("highly unreliable", conditionMessage(w)))
})

test_that("G = 2 boundary computes with df=1 and triggers G<10 strong warning", {
  df <- make_cluster_data_g(G = 2)
  fit <- lm(y ~ D + x, data = df)
  # Capture result while allowing warning to propagate
  warnings_caught <- list()
  res <- withCallingHandlers(
    compute_vce(fit, vce = "cluster", cluster = df$cl),
    lwdid_small_sample = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  # Verify warning was triggered with correct message

  expect_true(length(warnings_caught) >= 1L)
  expect_match(conditionMessage(warnings_caught[[1L]]), "highly unreliable")
  # Verify result
  expect_equal(res$df, 1L)
  expect_equal(res$n_clusters, 2L)
})

# ============================================================================
# Group 6: Imbalance and near-zero SE warnings
# ============================================================================

test_that("CV > 1.0 triggers lwdid_small_sample warning with 'unbalanced'", {
  # Use G=25 (>= 20) to avoid small-sample warning mixing
  # Highly unbalanced: 24 clusters of size 1, 1 cluster of size 100
  set.seed(42)
  n <- 24 * 1 + 100
  cluster <- c(rep(1:24, each = 1), rep(25L, 100))
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  expect_warning(
    compute_vce(fit, vce = "cluster", cluster = cluster),
    regexp = "unbalanced",
    class = "lwdid_small_sample"
  )
})

test_that("perfect fit triggers lwdid_numerical warning", {
  # y is exactly a linear function of D and x (no noise)
  set.seed(42)
  n <- 30
  d <- c(rep(0, 15), rep(1, 15))
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x  # no noise -> perfect fit -> SE ~ 0
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  # Expect lwdid_numerical warning (may also get R's "essentially perfect fit" warning)
  expect_warning(
    compute_vce(fit, vce = NULL),
    class = "lwdid_numerical"
  )
})

test_that("normal data with G >= 20 balanced clusters produces no warnings", {
  df <- make_cluster_data_g(G = 25)
  fit <- lm(y ~ D + x, data = df)
  expect_no_warning(
    compute_vce(fit, vce = "cluster", cluster = df$cl)
  )
})

test_that("model without 'D' coefficient safely skips near-zero SE check (AC-32)", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, x = x)
  fit <- lm(y ~ x, data = df)
  # Should not error or warn about near-zero SE for D
  res <- compute_vce(fit, vce = NULL)
  expect_true(is.matrix(res$vcov))
  expect_identical(res$vce_type, "homoskedastic")
})

# ============================================================================
# Group 7: verbose / WarningRegistry behavior
# ============================================================================

test_that("verbose parameter accepted by compute_vce() signature", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  # All three verbose values should be accepted without error
  expect_no_error(compute_vce(fit, vce = NULL, verbose = "quiet"))
  expect_no_error(compute_vce(fit, vce = NULL, verbose = "default"))
  expect_no_error(compute_vce(fit, vce = NULL, verbose = "verbose"))
})

test_that("warnings still fire with verbose='quiet'", {
  # Use G=5 to trigger small sample warning
  df <- make_cluster_data_g(G = 5)
  fit <- lm(y ~ D + x, data = df)
  # Even with verbose="quiet", warn_lwdid() still fires
  expect_warning(
    compute_vce(fit, vce = "cluster", cluster = df$cl, verbose = "quiet"),
    class = "lwdid_small_sample"
  )
})

# ============================================================================
# Group 8: Cluster variable type support
# ============================================================================

test_that("character cluster variable works correctly", {
  set.seed(42)
  n <- 100L
  G <- 20L
  cl_char <- rep(paste0("cl_", 1:G), each = n / G)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "cluster", cluster = cl_char)
  expect_true(is.matrix(res$vcov))
  expect_identical(res$vce_type, "cluster")
  expect_identical(res$n_clusters, G)
})

test_that("numeric cluster variable works correctly", {
  df <- make_cluster_data(n = 100, G = 20)
  fit <- lm(y ~ D + x, data = df)
  res <- compute_vce(fit, vce = "cluster", cluster = df$cl)
  expect_true(is.matrix(res$vcov))
  expect_identical(res$vce_type, "cluster")
  expect_identical(res$n_clusters, 20L)
})

test_that("character and numeric cluster variables produce same VCE result", {
  set.seed(42)
  n <- 100
  G <- 20
  cl_num <- rep(1:G, each = n / G)
  cl_char <- paste0("cl_", cl_num)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  res_num <- compute_vce(fit, vce = "cluster", cluster = cl_num)
  res_char <- compute_vce(fit, vce = "cluster", cluster = cl_char)
  expect_equal(res_num$vcov, res_char$vcov, tolerance = 1e-12)
})

# ============================================================================
# Group 9: HC hierarchy validation
# ============================================================================

test_that("HC1 > HC0 strictly: se_hc1 == se_hc0 * sqrt(n / (n - p)) exactly", {
  df <- make_test_data(n = 50)
  fit <- lm(y ~ D + x, data = df)
  n <- nobs(fit)
  p <- length(coef(fit))

  res_hc0 <- compute_vce(fit, vce = "hc0")
  res_hc1 <- compute_vce(fit, vce = "hc1")

  se_hc0 <- sqrt(diag(res_hc0$vcov))
  se_hc1 <- sqrt(diag(res_hc1$vcov))

  expected_se_hc1 <- se_hc0 * sqrt(n / (n - p))
  expect_equal(se_hc1, expected_se_hc1, tolerance = 1e-12)
})

test_that("small sample general pattern: HC3 >= HC2 >= HC1 >= HC0 (heteroskedastic data)", {
  # Construct heteroskedastic data
  set.seed(42)
  n <- 30
  d <- c(rep(0, 15), rep(1, 15))
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = abs(x) + 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  se_hc0 <- sqrt(diag(compute_vce(fit, vce = "hc0")$vcov))
  se_hc1 <- sqrt(diag(compute_vce(fit, vce = "hc1")$vcov))
  se_hc2 <- sqrt(diag(compute_vce(fit, vce = "hc2")$vcov))
  se_hc3 <- sqrt(diag(compute_vce(fit, vce = "hc3")$vcov))

  # HC1 >= HC0 (always true by construction)
  expect_true(all(se_hc1 >= se_hc0 - 1e-15))
  # HC2 >= HC1 (generally true for heteroskedastic data)
  expect_true(all(se_hc2 >= se_hc1 - 1e-15))
  # HC3 >= HC2 (generally true for heteroskedastic data)
  expect_true(all(se_hc3 >= se_hc2 - 1e-15))
})

# ============================================================================
# Group 10: compute_hc_vce() input validation and warnings
# ============================================================================

test_that("invalid HC type 'HC5' raises lwdid_invalid_parameter error", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  expect_error(
    compute_hc_vce(fit, type = "HC5"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("HC1/HC3 with extreme small sample (N_treated=1) triggers lwdid_small_sample warning", {
  set.seed(42)
  n <- 10
  d <- c(rep(0, 9), 1)  # N_treated = 1
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  expect_warning(
    compute_hc_vce(fit, type = "hc1"),
    class = "lwdid_small_sample"
  )
  expect_warning(
    compute_hc_vce(fit, type = "hc3"),
    class = "lwdid_small_sample"
  )
})

test_that("HC2-HC4 with extreme high leverage triggers lwdid_numerical warning", {
  set.seed(42)
  n <- 20
  x <- c(rep(0, n - 1), 1000)  # one extreme point
  d <- c(rep(0, n / 2), rep(1, n / 2))
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  # Verify the extreme leverage condition holds
  expect_true(max(hatvalues(fit)) > 0.99)
  # HC2, HC3, HC4 should all trigger the high leverage warning
  expect_warning(
    compute_hc_vce(fit, type = "hc2"),
    class = "lwdid_numerical"
  )
  expect_warning(
    compute_hc_vce(fit, type = "hc3"),
    class = "lwdid_numerical"
  )
  expect_warning(
    compute_hc_vce(fit, type = "hc4"),
    class = "lwdid_numerical"
  )
})

# ============================================================================
# Group 11: HC3-jackknife relationship
# ============================================================================

test_that("V_jack = ((N-1)/N) * V_HC3 exactly", {
  set.seed(42)
  n <- 10
  x <- rnorm(n)
  d <- c(rep(0, 5), rep(1, 5))
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  # Delete-one jackknife centered on beta_hat (full-sample estimate)
  beta_hat <- coef(fit)
  p <- length(beta_hat)
  beta_jack <- matrix(NA, n, p)
  for (i in 1:n) {
    fit_i <- lm(y ~ D + x, data = df[-i, ])
    beta_jack[i, ] <- coef(fit_i)
  }
  centered <- sweep(beta_jack, 2, beta_hat)
  V_jack <- ((n - 1) / n) * crossprod(centered)

  # HC3 VCE via sandwich
  V_hc3 <- sandwich::vcovHC(fit, type = "HC3")

  # Relationship: V_jack(centered on beta_hat) = ((N-1)/N) * V_HC3
  expect_equal(unname(V_jack), unname(((n - 1) / n) * V_hc3), tolerance = 1e-10)
})

# ============================================================================
# Group 12: Rank-deficient design matrix
# ============================================================================

test_that("rank-deficient design matrix does not crash compute_vce()", {
  set.seed(42)
  n <- 30
  x1 <- rnorm(n)
  x2 <- 2 * x1  # perfectly collinear
  d <- rbinom(n, 1, 0.5)
  y <- 1 + d + x1 + rnorm(n)
  df <- data.frame(y = y, D = d, x1 = x1, x2 = x2)
  fit <- lm(y ~ D + x1 + x2, data = df)
  # Should not crash
  res <- compute_vce(fit, vce = NULL)
  expect_true(is.matrix(res$vcov))
})

# ============================================================================
# Group 13: HC0-HC4 loop consistency against sandwich::vcovHC
# ============================================================================

test_that("HC0-HC4 all match sandwich::vcovHC exactly (< 1e-12)", {
  set.seed(42)
  n <- 30
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  
  for (hc_type in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    res <- compute_vce(fit, vce = hc_type)
    expected <- sandwich::vcovHC(fit, type = toupper(hc_type))
    expect_equal(res$vcov, expected, tolerance = 1e-12,
                 label = paste("compute_vce vcov for", hc_type),
                 expected.label = paste("sandwich::vcovHC for", toupper(hc_type)))
  }
})

# ============================================================================
# Group 14: Cluster VCE consistency against sandwich::vcovCL
# ============================================================================

test_that("cluster VCE matches sandwich::vcovCL(type='HC1') exactly (< 1e-12)", {
  set.seed(42)
  n <- 100L
  G <- 10L
  cluster <- rep(1:G, each = n / G)
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  
  # Suppress small sample warning (G=10 triggers informational warning)
  res <- suppressWarnings(compute_vce(fit, vce = "cluster", cluster = cluster))
  expected <- sandwich::vcovCL(fit, cluster = cluster, type = "HC1")
  
  expect_equal(res$vcov, expected, tolerance = 1e-12)
  expect_equal(res$df, G - 1L)
  expect_identical(res$n_clusters, G)
})

# ============================================================================
# Group 15: robust alias consistency
# ============================================================================

test_that("vce='robust' produces identical vcov to vce='hc1'", {
  set.seed(42)
  n <- 50
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  
  res_robust <- compute_vce(fit, vce = "robust")
  res_hc1 <- compute_vce(fit, vce = "hc1")
  
  expect_identical(res_robust$vcov, res_hc1$vcov)
  expect_identical(res_robust$vce_type, res_hc1$vce_type)
})


# ============================================================================
# Group 16: compute_hc_vce() sandwich consistency (E3-02.2)
# ============================================================================

test_that("compute_hc_vce() HC0-HC4 SE matches sandwich::vcovHC exactly (< 1e-12)", {
  set.seed(42)
  n <- 30
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = abs(x) + 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  for (hc in paste0("HC", 0:4)) {
    res <- compute_hc_vce(fit, type = hc)
    expected_vcov <- sandwich::vcovHC(fit, type = hc)

    # vcov matches sandwich
    expect_equal(
      res$vcov, expected_vcov,
      tolerance = 1e-12,
      label = paste("compute_hc_vce vcov for", hc),
      expected.label = paste("sandwich::vcovHC for", hc)
    )

    # se matches sandwich
    expected_se <- sqrt(diag(expected_vcov))
    expect_equal(
      res$se, expected_se,
      tolerance = 1e-12,
      label = paste("compute_hc_vce se for", hc),
      expected.label = paste("sqrt(diag(sandwich)) for", hc)
    )

    # se == sqrt(diag(vcov)) exactly
    expect_equal(
      res$se, sqrt(diag(res$vcov)),
      tolerance = 1e-15,
      label = paste("se vs sqrt(diag(vcov)) for", hc)
    )
  }
})

test_that("compute_hc_vce() robust alias via dispatcher matches HC1 exactly", {
  set.seed(42)
  n <- 30
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = abs(x) + 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  res_robust <- compute_vce(fit, vce = "robust")
  res_hc1 <- compute_hc_vce(fit, type = "HC1")

  expect_identical(res_robust$vcov, res_hc1$vcov)
})

test_that("compute_hc_vce() returns list with exactly vcov and se fields", {
  set.seed(42)
  n <- 30
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = abs(x) + 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  res <- compute_hc_vce(fit, type = "HC3")

  # Exactly two fields: vcov and se
  expect_identical(sort(names(res)), c("se", "vcov"))

  # vcov is a matrix
  expect_true(is.matrix(res$vcov))

  # se is a numeric vector
  expect_true(is.numeric(res$se))
  expect_true(is.vector(res$se))

  # Dimensions are consistent
  expect_equal(length(res$se), ncol(res$vcov))
  expect_equal(length(res$se), nrow(res$vcov))
})


# ============================================================================
# Group 17: HC hierarchy and large sample convergence (E3-02.3)
# ============================================================================

test_that("HC1 > HC0 strict scalar relationship via compute_hc_vce() directly", {
  set.seed(123)
  n <- 50
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = abs(x) + 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  p <- length(coef(fit))

  res_hc0 <- compute_hc_vce(fit, type = "HC0")
  res_hc1 <- compute_hc_vce(fit, type = "HC1")

  # SE relationship: se_hc1 == se_hc0 * sqrt(n / (n - p))
  expected_se <- res_hc0$se * sqrt(n / (n - p))
  expect_equal(res_hc1$se, expected_se, tolerance = 1e-12)

  # Vcov relationship: vcov_hc1 == (n / (n - p)) * vcov_hc0
  expected_vcov <- (n / (n - p)) * res_hc0$vcov
  expect_equal(res_hc1$vcov, expected_vcov, tolerance = 1e-12)
})

test_that("small sample HC hierarchy via compute_hc_vce(): HC3 >= HC2 >= HC1 >= HC0", {
  set.seed(123)
  n <- 20
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  # Strongly heteroskedastic data
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = (abs(x) + 0.1)^2)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  se_hc0 <- compute_hc_vce(fit, type = "HC0")$se
  se_hc1 <- compute_hc_vce(fit, type = "HC1")$se
  se_hc2 <- compute_hc_vce(fit, type = "HC2")$se
  se_hc3 <- compute_hc_vce(fit, type = "HC3")$se

  # Verify HC3 >= HC2 >= HC1 >= HC0 for ALL coefficients
  expect_true(all(se_hc1 >= se_hc0 - 1e-10))
  expect_true(all(se_hc2 >= se_hc1 - 1e-10))
  expect_true(all(se_hc3 >= se_hc2 - 1e-10))
})

test_that("large sample convergence: HC0-HC4 SE relative diff < 5% at N=1000", {
  set.seed(42)
  n <- 1000
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = abs(x) + 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)
  p <- length(coef(fit))

  # Compute SE for HC0-HC4 via compute_hc_vce() directly
  se_list <- list()
  for (hc in paste0("HC", 0:4)) {
    se_list[[hc]] <- compute_hc_vce(fit, type = hc)$se
  }

  # For each pair of HC types, verify relative SE difference < 5%
  hc_names <- paste0("HC", 0:4)
  for (i in seq_along(hc_names)) {
    for (j in seq_along(hc_names)) {
      if (i < j) {
        se_a <- se_list[[hc_names[i]]]
        se_b <- se_list[[hc_names[j]]]
        rel_diff <- abs(se_a - se_b) / se_a
        expect_true(
          all(rel_diff < 0.05),
          label = paste(
            hc_names[i], "vs", hc_names[j],
            "max rel diff =", round(max(rel_diff), 4)
          )
        )
      }
    }
  }

  # Verify max(hatvalues) is small (< 0.05), confirming h_ii = O(1/N)
  h <- hatvalues(fit)
  expect_true(max(h) < 0.05)

  # Verify mean(hatvalues) ≈ p/n
  expect_equal(mean(h), p / n, tolerance = 1e-10)
})

# ============================================================================
# Group 18: Mathematical properties verification (E3-02.4)
# ============================================================================

test_that("HC2 unbiasedness: homoskedastic data HC2 SE close to OLS SE (< 20%)", {
  set.seed(42)
  n <- 50
  d <- c(rep(0, 25), rep(1, 25))
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = 1)  # known sigma^2 = 1
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  se_ols <- sqrt(diag(vcov(fit)))
  se_hc2 <- compute_hc_vce(fit, type = "HC2")$se

  # Under homoskedasticity, E[e_i^2/(1-h_ii)] = sigma^2 (REQ-3215)
  # so HC2 SE should be close to OLS SE
  rel_diff <- abs(se_hc2 - se_ols) / se_ols
  expect_true(
    all(rel_diff < 0.20),
    label = paste("HC2 vs OLS relative diff: max =", round(max(rel_diff), 4))
  )
})

test_that("HC3-jackknife relationship via compute_hc_vce(): V_jack = ((N-1)/N) * V_HC3", {
  set.seed(99)
  n <- 15
  d <- c(rep(0, 8), rep(1, 7))
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = 0.8)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  # HC3 via compute_hc_vce()
  V_hc3 <- compute_hc_vce(fit, type = "HC3")$vcov

  # Manual delete-one jackknife
  beta_hat <- coef(fit)
  p <- length(beta_hat)
  beta_jack <- matrix(NA, n, p)
  for (i in 1:n) {
    fit_i <- lm(y ~ D + x, data = df[-i, ])
    beta_jack[i, ] <- coef(fit_i)
  }
  centered <- sweep(beta_jack, 2, beta_hat)
  V_jack <- ((n - 1) / n) * crossprod(centered)

  # Relationship: V_jack(centered on beta_hat) = ((N-1)/N) * V_HC3
  expect_equal(
    unname(V_jack),
    unname(((n - 1) / n) * V_hc3),
    tolerance = 1e-10
  )
})

test_that("Sherman-Morrison formula: beta_{(-i)} = beta - (X'X)^{-1} x_i e_i / (1 - h_ii)", {
  set.seed(42)
  n <- 10
  d <- c(rep(0, 5), rep(1, 5))
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  beta_hat <- coef(fit)
  X <- model.matrix(fit)
  e <- resid(fit)
  h <- hatvalues(fit)
  XtX_inv <- solve(crossprod(X))

  for (i in 1:n) {
    # Actual delete-one estimate
    beta_loo <- coef(lm(y ~ D + x, data = df[-i, ]))

    # Sherman-Morrison prediction
    beta_sm <- beta_hat - as.vector(XtX_inv %*% X[i, ] * e[i] / (1 - h[i]))

    expect_equal(
      unname(beta_loo), unname(beta_sm),
      tolerance = 1e-10,
      label = paste("Sherman-Morrison for observation", i)
    )
  }
})

test_that("HC4 d_i behavior: high leverage d_i > 3.5, normal d_i < 1.5", {
  set.seed(42)
  n <- 20
  x <- c(rnorm(n - 1, mean = 0, sd = 1), 100)  # last point is extreme
  d <- c(rep(0, n / 2), rep(1, n / 2))
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  h <- hatvalues(fit)
  h_bar <- mean(h)
  d_i <- pmin(4, h / h_bar)

  # High-leverage observation (last one) should have d_i > 3.5
  expect_true(d_i[n] > 3.5,
              label = paste("High-leverage d_i =", round(d_i[n], 4)))

  # Normal observations should have d_i < 1.5
  expect_true(all(d_i[-n] < 1.5),
              label = paste("Max normal d_i =", round(max(d_i[-n]), 4)))

  # Natural non-negativity: d_i >= 0 for all
  expect_true(all(d_i >= 0))
})

test_that("Leverage values: sum(h_ii) = p and all h_ii in [0, 1)", {
  set.seed(42)
  n <- 50
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  h <- hatvalues(fit)
  p <- length(coef(fit))

  # sum(h_ii) = p (trace of hat matrix equals number of parameters)
  expect_equal(sum(h), p, tolerance = 1e-10)

  # All h_ii >= 0 (non-negative)
  expect_true(all(h >= 0))

  # All h_ii < 1 (strict, since N > p)
  expect_true(all(h < 1))

  # Lower bound: h_ii >= 1/n for models with intercept
  expect_true(all(h >= 1 / n - 1e-10))
})

# ============================================================================
# Group 19: compute_hc_vce() input validation and warnings (E3-02.5)
# ============================================================================

test_that("invalid HC type 'invalid' raises lwdid_invalid_parameter via compute_hc_vce()", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  expect_error(
    compute_hc_vce(fit, type = "invalid"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("empty string type raises lwdid_invalid_parameter via compute_hc_vce()", {
  df <- make_test_data()
  fit <- lm(y ~ D + x, data = df)
  expect_error(
    compute_hc_vce(fit, type = ""),
    class = "lwdid_invalid_parameter"
  )
})

test_that("HC0 does NOT trigger small sample warning even with N_treated=1", {
  set.seed(42)
  n <- 10
  d <- c(rep(0, 9), 1)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  # HC0 should NOT trigger lwdid_small_sample (only HC1/HC3 do)
  warnings_caught <- list()
  res <- withCallingHandlers(
    compute_hc_vce(fit, type = "HC0"),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  small_sample_warnings <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_caught
  )
  expect_length(small_sample_warnings, 0L)
  expect_true(is.matrix(res$vcov))
})

test_that("HC2 does NOT trigger small sample warning even with N_treated=1", {
  set.seed(42)
  n <- 10
  d <- c(rep(0, 9), 1)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  # HC2 may trigger lwdid_numerical (high leverage), but NOT lwdid_small_sample
  warnings_caught <- list()
  res <- withCallingHandlers(
    compute_hc_vce(fit, type = "HC2"),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  small_sample_warnings <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_caught
  )
  expect_length(small_sample_warnings, 0L)
  expect_true(is.matrix(res$vcov))
})

test_that("HC0 and HC1 do NOT trigger high leverage warning even with max(h_ii) > 0.99", {
  set.seed(42)
  n <- 20
  x <- c(rep(0, n - 1), 1000)
  d <- c(rep(0, n / 2), rep(1, n / 2))
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  # Confirm extreme leverage condition

  expect_true(max(hatvalues(fit)) > 0.99)

  # HC0: no lwdid_numerical warning expected
  warnings_hc0 <- list()
  withCallingHandlers(
    compute_hc_vce(fit, type = "HC0"),
    warning = function(w) {
      warnings_hc0[[length(warnings_hc0) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  numerical_hc0 <- Filter(
    function(w) inherits(w, "lwdid_numerical"),
    warnings_hc0
  )
  expect_length(numerical_hc0, 0L)

  # HC1: no lwdid_numerical warning expected
  warnings_hc1 <- list()
  withCallingHandlers(
    compute_hc_vce(fit, type = "HC1"),
    warning = function(w) {
      warnings_hc1[[length(warnings_hc1) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  numerical_hc1 <- Filter(
    function(w) inherits(w, "lwdid_numerical"),
    warnings_hc1
  )
  expect_length(numerical_hc1, 0L)
})

test_that("model without D variable safely skips small sample check in compute_hc_vce()", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, x = x)
  fit <- lm(y ~ x, data = df)

  # HC1 should not error or warn about small sample (no D variable)
  warnings_hc1 <- list()
  res_hc1 <- withCallingHandlers(
    compute_hc_vce(fit, type = "HC1"),
    warning = function(w) {
      warnings_hc1[[length(warnings_hc1) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  small_sample_hc1 <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_hc1
  )
  expect_length(small_sample_hc1, 0L)
  expect_true(is.matrix(res_hc1$vcov))

  # HC3 should not error or warn about small sample (no D variable)
  warnings_hc3 <- list()
  res_hc3 <- withCallingHandlers(
    compute_hc_vce(fit, type = "HC3"),
    warning = function(w) {
      warnings_hc3[[length(warnings_hc3) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  small_sample_hc3 <- Filter(
    function(w) inherits(w, "lwdid_small_sample"),
    warnings_hc3
  )
  expect_length(small_sample_hc3, 0L)
  expect_true(is.matrix(res_hc3$vcov))
})

test_that("default type parameter is HC3", {
  set.seed(42)
  n <- 30
  d <- rbinom(n, 1, 0.5)
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  res_default <- compute_hc_vce(fit)
  res_hc3 <- compute_hc_vce(fit, type = "HC3")

  expect_identical(res_default$vcov, res_hc3$vcov)
  expect_identical(res_default$se, res_hc3$se)
})

# ============================================================================
# Group 20: Smoking dataset benchmark (E3-02.6)
# ============================================================================
# Synthetic panel data inspired by the California Proposition 99
# anti-smoking study (Abadie, Diamond & Hainmueller, 2010).
# 39 states × 20 periods = 780 obs, two-way FE (state + period).
# California (state 1) treated after period 10.
# ============================================================================

#' Generate smoking-style DiD panel data
make_smoking_data <- function(seed = 2010) {
  set.seed(seed)
  n_states  <- 39L
  n_periods <- 20L
  n <- n_states * n_periods  # 780

  state  <- rep(1:n_states, each = n_periods)
  period <- rep(1:n_periods, times = n_states)

  # California is state 1, treated after period 10 (1980)
  D <- as.integer(state == 1L & period > 10L)

  # State-level heterogeneity
  state_effect <- rep(rnorm(n_states, sd = 5), each = n_periods)

  # Time trend
  time_trend <- rep(seq(0, 2, length.out = n_periods), times = n_states)

  # Heteroskedastic errors (larger for California)
  sigma <- ifelse(state == 1L, 3, 1)
  y <- 100 + state_effect + 2 * time_trend - 15 * D + rnorm(n, sd = sigma)

  data.frame(
    y      = y,
    D      = D,
    state  = factor(state),
    period = factor(period)
  )
}

test_that("Smoking-style DiD data: HC3 SE matches sandwich exactly", {
  df  <- make_smoking_data(seed = 2010)
  fit <- lm(y ~ D + state + period, data = df)

  # HC3 via compute_hc_vce()
  res <- suppressWarnings(compute_hc_vce(fit, type = "HC3"))

  # HC3 via sandwich directly
  expected_vcov <- sandwich::vcovHC(fit, type = "HC3")

  # Full vcov matrix match

  expect_equal(res$vcov, expected_vcov, tolerance = 1e-12)

  # SE vector match
  expected_se <- sqrt(diag(expected_vcov))
  expect_equal(res$se, expected_se, tolerance = 1e-12)

  # Extract SE for the treatment coefficient "D"
  d_idx <- which(names(coef(fit)) == "D")
  se_D  <- res$se[d_idx]

  # SE for D should be positive and in a reasonable range
  expect_true(se_D > 0.1,
              label = paste("SE(D) =", round(se_D, 6), "> 0.1"))
  expect_true(se_D < 50,
              label = paste("SE(D) =", round(se_D, 6), "< 50"))
})

test_that("Smoking-style data: HC hierarchy holds", {
  df  <- make_smoking_data(seed = 2010)
  fit <- lm(y ~ D + state + period, data = df)

  # Compute HC0-HC3 SE for the "D" coefficient
  d_idx <- which(names(coef(fit)) == "D")

  se_hc0 <- suppressWarnings(compute_hc_vce(fit, type = "HC0"))$se[d_idx]
  se_hc1 <- suppressWarnings(compute_hc_vce(fit, type = "HC1"))$se[d_idx]
  se_hc2 <- suppressWarnings(compute_hc_vce(fit, type = "HC2"))$se[d_idx]
  se_hc3 <- suppressWarnings(compute_hc_vce(fit, type = "HC3"))$se[d_idx]

  # HC1 > HC0 strictly (finite sample correction n/(n-p))
  expect_true(se_hc1 > se_hc0,
              label = paste("HC1", round(se_hc1, 8), "> HC0", round(se_hc0, 8)))

  # HC3 >= HC2 >= HC1 (with tolerance for numerical noise)
  expect_true(se_hc2 >= se_hc1 - 1e-10,
              label = paste("HC2", round(se_hc2, 8), ">= HC1", round(se_hc1, 8)))
  expect_true(se_hc3 >= se_hc2 - 1e-10,
              label = paste("HC3", round(se_hc3, 8), ">= HC2", round(se_hc2, 8)))
})

test_that("Smoking-style data: leverage properties", {
  df  <- make_smoking_data(seed = 2010)
  fit <- lm(y ~ D + state + period, data = df)

  h <- hatvalues(fit)
  p <- length(coef(fit))

  # sum(h_ii) = p (trace of hat matrix equals number of parameters)
  # Larger tolerance because p is large with many FE dummies (~57)
  expect_equal(sum(h), p, tolerance = 1e-8,
               label = paste("sum(h) =", round(sum(h), 6), "vs p =", p))

  # All h_ii >= 0
  expect_true(all(h >= 0),
              label = "All leverage values non-negative")

  # All h_ii < 1 (strict, since N >> p)
  expect_true(all(h < 1),
              label = paste("max(h) =", round(max(h), 6), "< 1"))

  # Stata benchmark values (to be filled via Stata MCP): HC3 SE for D = ???
})
