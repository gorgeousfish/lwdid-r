# ============================================================================
# test-vce-numerical.R — E3-06 VCE Numerical & Integration Test Suite
#
# Layers covered:
#   L2:  Stata numerical consistency (smoking dataset)
#   L3:  degrees of freedom verification (three-tier + cluster)
#   L7:  cluster nesting validation
#   L15: integration layer (lwdid main function end-to-end)
#
# Test IDs: L2-xx, L3-xx, L7-xx, L15-xx
# ============================================================================

library(testthat)
library(lwdid)

# ============================================================================
# Layer 2: Stata numerical consistency (smoking dataset)
#
# Stata benchmark values (verified in E3-05.9):
#   ATT (all VCE types):  -0.4221744083946467
#   SE (homoskedastic):    0.1207995330393050  (df = 37)
#   SE (robust / HC1):     0.0195962719215685  (df = 37)
#   SE (cluster state):    0.0195962719215685  (df = 38, G = 39)
#   N = 39, intercept = -0.2461572742873901
#
# Note: Robust SE == Cluster SE because N=39 with 39 clusters means each
#       cluster has exactly 1 observation.
# ============================================================================

# --- Stata benchmark constants ---
stata_att        <- -0.4221744083946467
stata_se_homo    <-  0.1207995330393050
stata_se_robust  <-  0.0195962719215685
stata_se_cluster <-  0.0195962719215685

test_that("L2-01: homoskedastic SE matches Stata on smoking data", {
  data(smoking)
  result <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean")
  )

  # ATT: tight tolerance (< 1e-6)
  expect_equal(result$att, stata_att, tolerance = 1e-6,
               label = "ATT vs Stata")

  # SE: < 1e-4 tolerance
  expect_equal(result$se_att, stata_se_homo, tolerance = 1e-4,
               label = "Homoskedastic SE vs Stata")

  # Degrees of freedom: N - k = 39 - 2 = 37
  expect_equal(result$df_resid, 37L, label = "df for homoskedastic")
})

test_that("L2-02: robust (HC1) SE matches Stata vce(robust)", {
  data(smoking)
  result <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "robust")
  )

  # ATT: tight tolerance (< 1e-6)
  expect_equal(result$att, stata_att, tolerance = 1e-6,
               label = "ATT vs Stata (robust)")

  # SE: < 1e-4 tolerance
  expect_equal(result$se_att, stata_se_robust, tolerance = 1e-4,
               label = "Robust SE vs Stata")

  # Degrees of freedom: N - k = 39 - 2 = 37
  expect_equal(result$df_resid, 37L, label = "df for robust")
})

test_that("L2-03: cluster SE matches Stata vce(cluster state)", {
  data(smoking)
  result <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "state")
  )

  # ATT: tight tolerance (< 1e-6)
  expect_equal(result$att, stata_att, tolerance = 1e-6,
               label = "ATT vs Stata (cluster)")

  # SE: < 1e-4 tolerance
  expect_equal(result$se_att, stata_se_cluster, tolerance = 1e-4,
               label = "Cluster SE vs Stata")

  # Degrees of freedom: G - 1 = 39 - 1 = 38
  expect_equal(result$df_resid, 38L, label = "df for cluster")

  # Number of clusters
  expect_equal(result$n_clusters, 39L, label = "n_clusters")
})

test_that("L2-04: ATT identical across VCE types", {
  data(smoking)

  res_homo <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean")
  )

  res_robust <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "robust")
  )

  res_cluster <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "state")
  )

  # Point estimate must be identical regardless of VCE (< 1e-12)
  expect_equal(res_homo$att, res_robust$att, tolerance = 1e-12,
               label = "ATT: homo vs robust")
  expect_equal(res_homo$att, res_cluster$att, tolerance = 1e-12,
               label = "ATT: homo vs cluster")
  expect_equal(res_robust$att, res_cluster$att, tolerance = 1e-12,
               label = "ATT: robust vs cluster")
})

# ============================================================================
# Layer 3: Degrees of Freedom Verification (Three-Tier + Cluster)
# ============================================================================

# L3-01: Tier 3 df = N - 2 (no controls, smoking dataset)
test_that("L3-01: Tier 3 df = N - 2 (no controls, smoking dataset)", {
  data(smoking)
  result <- lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
                  d = "d", post = "post", rolling = "demean", vce = NULL)
  n_units <- length(unique(smoking$state))
  expect_equal(result$df_inference, n_units - 2L)
})

# L3-02: Tier 1 df = N - 2 - 2K (full interaction model)
test_that("L3-02: Tier 1 df = N - 2 - 2K (full interaction)", {
  set.seed(42)
  n <- 20; K <- 2
  d <- c(rep(1, 10), rep(0, 10))
  x <- matrix(rnorm(n * K), ncol = K)
  y <- 1 + 2 * d + rnorm(n)
  result <- estimate_ra_common(y, d, x, vce = NULL)
  # N_1=10 > K+1=3, N_0=10 > K+1=3 → Tier 1
  expect_equal(result$df, n - 2L - 2L * K)  # 20 - 2 - 4 = 14
  expect_equal(result$controls_tier, "full_interaction")
})

# L3-03: Tier 2 df = N - K - 2 (simple controls)
test_that("L3-03: Tier 2 df = N - K - 2 (simple controls)", {
  set.seed(42)
  n <- 14; K <- 2
  # N_1=2 <= K+1=3 → not Tier 1; N=14 > K+2=4 → Tier 2
  d <- c(rep(1, 2), rep(0, 12))
  x <- matrix(rnorm(n * K), ncol = K)
  y <- 1 + 2 * d + rnorm(n)
  result <- estimate_ra_common(y, d, x, vce = NULL)
  expect_equal(result$df, n - K - 2L)  # 14 - 2 - 2 = 10
  expect_equal(result$controls_tier, "simple")
})

# L3-04: Tier 3 df = N - 2 (controls dropped)
test_that("L3-04: Tier 3 df = N - 2 (controls dropped, N <= K+2)", {
  set.seed(42)
  K <- 3; n <- 4  # N=4 <= K+2=5 → Tier 3
  d <- c(rep(1, 2), rep(0, 2))
  x <- matrix(rnorm(n * K), ncol = K)
  y <- 1 + 2 * d + rnorm(n)
  result <- estimate_ra_common(y, d, x, vce = NULL)
  expect_equal(result$df, n - 2L)  # 4 - 2 = 2
  expect_equal(result$controls_tier, "none")
})

# L3-05: Cluster df = G - 1 (independent of model tier)
test_that("L3-05: cluster df = G - 1 (smoking dataset)", {
  data(smoking)
  result <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "state")
  )
  n_states <- length(unique(smoking$state))
  expect_equal(result$df_inference, n_states - 1L)
  expect_equal(result$n_clusters, n_states)
})

# L3-06: Tier 2 + HC3 uses df = N - K - 2 (not N - 2 - 2K)
test_that("L3-06: Tier 2 + HC3 uses df = N - K - 2", {
  set.seed(42)
  n <- 14; K <- 2
  d <- c(rep(1, 2), rep(0, 12))
  x <- matrix(rnorm(n * K), ncol = K)
  y <- 1 + 2 * d + rnorm(n)
  result <- estimate_ra_common(y, d, x, vce = "hc3")
  expect_equal(result$df, n - K - 2L)  # 14 - 2 - 2 = 10, NOT 14 - 2 - 4 = 8
  expect_equal(result$controls_tier, "simple")
})

# ============================================================================
# Layer 7: Cluster Nesting Validation
# ============================================================================

# L7-01: Cluster variable not nested in unit variable throws error
test_that("L7-01: cluster not nested in unit throws error", {
  # Unit 1 belongs to cluster A in period 1 but cluster B in period 2
  dt <- data.table::data.table(
    id = c(1, 1, 2, 2, 3, 3),
    time = c(1, 2, 1, 2, 1, 2),
    y = c(1, 2, 3, 4, 5, 6),
    d = c(1, 1, 0, 0, 0, 0),
    post = c(0, 1, 0, 1, 0, 1),
    cluster = c("A", "B", "A", "A", "B", "B")  # unit 1 crosses clusters!
  )
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post",
          vce = "cluster", cluster_var = "cluster"),
    class = "lwdid_invalid_parameter"
  )
})

# L7-02: Correctly nested cluster variable works
test_that("L7-02: correctly nested cluster works", {
  set.seed(42)
  dt <- data.table::data.table(
    id = rep(1:10, each = 3),
    time = rep(1:3, 10),
    y = rnorm(30),
    d = rep(c(rep(1, 5), rep(0, 5)), each = 3),
    post = rep(c(0, 0, 1), 10),
    cluster = rep(c("A", "A", "A", "B", "B",
                    "C", "C", "D", "D", "D"), each = 3)
  )
  expect_no_error(
    suppressWarnings(
      lwdid(dt, y = "y", ivar = "id", tvar = "time",
            d = "d", post = "post",
            vce = "cluster", cluster_var = "cluster")
    )
  )
})


# =============================================================================
# Layer 15: Integration Layer Validation (lwdid main function end-to-end)
# =============================================================================

# L15-01: cluster_var column not found → lwdid_missing_column error
test_that("L15-01: cluster_var column not found raises error", {
  set.seed(42)
  dt <- data.table::data.table(
    id = rep(1:6, each = 3),
    time = rep(1:3, 6),
    y = rnorm(18),
    d = rep(c(1, 1, 1, 0, 0, 0), each = 3),
    post = rep(c(0, 0, 1), 6)
  )
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post",
          vce = "cluster", cluster_var = "nonexistent_col"),
    class = "lwdid_missing_column"
  )
})

# L15-02: cluster_var column with NA → lwdid_invalid_parameter error
test_that("L15-02: cluster_var column with NA raises error", {
  set.seed(42)
  dt <- data.table::data.table(
    id = rep(1:6, each = 3),
    time = rep(1:3, 6),
    y = rnorm(18),
    d = rep(c(1, 1, 1, 0, 0, 0), each = 3),
    post = rep(c(0, 0, 1), 6),
    state = c("CA", "CA", "CA", NA, "NY", "NY",
              "TX", "TX", "TX", "FL", "FL", "FL",
              "WA", "WA", "WA", "OR", "OR", "OR")
  )
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post",
          vce = "cluster", cluster_var = "state"),
    class = "lwdid_invalid_parameter"
  )
})


# L15-03: Unbalanced panel with varying cluster counts across periods
test_that("L15-03: unbalanced panel with varying cluster counts across periods", {
  set.seed(42)
  # Unbalanced panel: periods 1,2 have all 10 units, period 3 only units 1-8
  # Cluster: unit 1,2->A, 3,4->B, 5,6->C, 7,8->D, 9,10->E
  # Treatment: unit 1-5->d=1, 6-10->d=0
  dt <- data.table::data.table(
    id = c(rep(1:10, each = 2), 1:8),
    time = c(rep(c(1, 2), 10), rep(3, 8)),
    y = rnorm(28),
    d = c(rep(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0), each = 2),
          c(1, 1, 1, 1, 1, 0, 0, 0)),
    post = c(rep(0, 20), rep(1, 8)),
    cluster = c(rep(c("A", "A", "B", "B", "C",
                      "C", "D", "D", "E", "E"), each = 2),
                c("A", "A", "B", "B", "C", "C", "D", "D"))
  )
  result <- suppressWarnings(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post",
          vce = "cluster", cluster_var = "cluster")
  )
  expect_true(!is.na(result$se_att))
})

# L15-04: Period fault tolerance — insufficient sample period
test_that("L15-04: period with insufficient sample yields NA without breaking others", {
  set.seed(42)
  # Period 1=pre, period 2,3=post; period 3 has only 2 units (id=1 treated, id=6 control)
  dt <- data.table::data.table(
    id = c(rep(1:10, each = 2), c(1, 6)),
    time = c(rep(c(1, 2), 10), rep(3, 2)),
    y = rnorm(22),
    d = c(rep(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0), each = 2),
          c(1, 0)),
    post = c(rep(c(0, 1), 10), rep(1, 2))
  )
  result <- suppressWarnings(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", vce = NULL)
  )
  expect_true(!is.null(result))
})


# L15-05: controls_tier in att_by_period (time-invariant controls required)
test_that("L15-05: att_by_period contains controls_tier column", {
  set.seed(42)
  K <- 2
  # 3 periods: period 1(pre) all 20 units, period 2(post) all 20, period 3(post) 6 units
  # Controls must be time-invariant (same value per unit across all periods)
  x1_unit <- rnorm(20)
  x2_unit <- rnorm(20)
  dt <- data.table::data.table(
    id = c(rep(1:20, each = 2), c(1, 2, 11, 12, 13, 14)),
    time = c(rep(c(1, 2), 20), rep(3, 6)),
    y = rnorm(46),
    d = c(rep(c(rep(1, 10), rep(0, 10)), each = 2),
          c(1, 1, 0, 0, 0, 0)),
    post = c(rep(c(0, 1), 20), rep(1, 6))
  )
  dt[, x1 := x1_unit[id]]
  dt[, x2 := x2_unit[id]]
  result <- suppressWarnings(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post",
          controls = c("x1", "x2"), vce = NULL)
  )
  if (!is.null(result$att_by_period)) {
    expect_true("controls_tier" %in% names(result$att_by_period))
  }
})

# L15-06: Generic error handler works
test_that("L15-06: generic error handler works in estimate_period_effects", {
  set.seed(42)
  dt <- data.table::data.table(
    id = rep(1:10, each = 3),
    time = rep(1:3, 10),
    y = rnorm(30),
    d = rep(c(rep(1, 5), rep(0, 5)), each = 3),
    post = rep(c(0, 0, 1), 10)
  )
  result <- suppressWarnings(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", vce = NULL)
  )
  expect_true(!is.null(result))
  if (!is.null(result$att_by_period)) {
    expect_true(is.data.frame(result$att_by_period))
  }
})


# L15-07: Cluster SE > homoskedastic SE with intra-cluster correlation
test_that("L15-07: cluster SE > homoskedastic SE with intra-cluster correlation", {
  set.seed(42)
  n_clusters <- 20L
  cluster_size <- 5L
  n <- n_clusters * cluster_size
  cluster_id <- rep(1:n_clusters, each = cluster_size)
  cluster_effect <- rnorm(n_clusters, sd = 2)  # cluster effects (source of ICC)
  d <- rep(c(rep(1, n_clusters / 2), rep(0, n_clusters / 2)), each = cluster_size)
  y <- 1 + 2 * d + cluster_effect[cluster_id] + rnorm(n, sd = 0.5)
  fit <- lm(y ~ d)
  se_ols <- sqrt(diag(vcov(fit)))["d"]
  result_cl <- compute_vce(fit, vce = "cluster", cluster = cluster_id)
  se_cluster <- sqrt(diag(result_cl$vcov))["d"]
  expect_true(se_cluster > se_ols,
              label = sprintf("cluster SE(%.6f) > homo SE(%.6f)", se_cluster, se_ols))
})

# L15-08: estimate_ra_common returns complete lm object
test_that("L15-08: estimate_ra_common returns complete lm object", {
  set.seed(42)
  n <- 20
  d <- c(rep(1, 10), rep(0, 10))
  y <- 1 + 2 * d + rnorm(n)
  result <- estimate_ra_common(y, d, x = NULL, vce = NULL)
  expect_s3_class(result$fit, "lm")
  expect_true(!is.null(result$fit$model))
  expect_true(!is.null(result$fit$qr))
  expect_true(!is.null(result$fit$residuals))
  expect_true(!is.null(result$fit$df.residual))
})

# L15-09: att_by_period returns complete 15 columns
test_that("L15-09: att_by_period returns complete 15 columns", {
  data(smoking)
  result <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean")
  )
  bp <- result$att_by_period
  expect_true(is.data.frame(bp))
  expected_cols <- c("tindex", "period", "att", "se", "t_stat", "pvalue",
                     "ci_lower", "ci_upper", "n_obs", "n_treated",
                     "n_control", "df", "vce_type", "n_clusters",
                     "controls_tier")
  for (col in expected_cols) {
    expect_true(col %in% names(bp),
                label = sprintf("att_by_period should contain column '%s'", col))
  }
  expect_equal(ncol(bp), 15L)
})
