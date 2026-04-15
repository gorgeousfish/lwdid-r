# ============================================================================
# test-vce-stata-consistency.R — E3-05.9 Stata Numerical Consistency Tests
#
# Verifies that lwdid R package VCE results match Stata benchmarks on the
# smoking dataset (39 US states, 31 years 1970-2000, California treated
# from 1989). With rolling="demean", the cross-sectional regression has
# N=39 observations.
#
# Stata benchmark values:
#   Regression: y_demean ~ d (Tier 3, no controls)
#   Homoskedastic OLS: ATT = -0.4221743, SE = 0.1207995, df = 37
#   Robust (HC1):      ATT = -0.4221743, SE = 0.0195963, df = 37
#   Cluster by state:  ATT = -0.4221743, SE = 0.0195963, df = 38, G = 39
#
# Test IDs: T5f-01 through T5f-04
# ============================================================================

# T5f-01: Homoskedastic SE matches Stata
test_that("T5f-01: homoskedastic ATT/SE/df match Stata benchmarks", {
  data(smoking)
  result <- lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
                  d = "d", post = "post", rolling = "demean")

  expect_equal(result$att, -0.4221743, tolerance = 1e-6)
  expect_equal(result$se_att, 0.1207995, tolerance = 1e-4)
  expect_equal(result$df_resid, 37L)
  expect_equal(result$nobs, 39L)
})

# T5f-02: Robust (HC1) SE matches Stata
test_that("T5f-02: robust HC1 SE matches Stata benchmark", {
  data(smoking)
  result <- lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
                  d = "d", post = "post", rolling = "demean",
                  vce = "robust")

  expect_equal(result$att, -0.4221743, tolerance = 1e-6)
  expect_equal(result$se_att, 0.0195963, tolerance = 1e-4)
  expect_equal(result$df_resid, 37L)
})

# T5f-03: Cluster SE matches Stata
test_that("T5f-03: cluster SE matches Stata benchmark", {
  data(smoking)
  result <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "state")
  )

  expect_equal(result$att, -0.4221743, tolerance = 1e-6)
  expect_equal(result$se_att, 0.0195963, tolerance = 1e-4)
  expect_equal(result$df_resid, 38L)
  expect_equal(result$n_clusters, 39L)
})

# T5f-04: Robust and cluster SEs identical (1 obs per cluster)
test_that("T5f-04: robust and cluster SEs identical for cross-sectional data", {
  data(smoking)
  res_robust <- lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
                      d = "d", post = "post", rolling = "demean",
                      vce = "robust")
  res_cluster <- suppressWarnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          vce = "cluster", cluster_var = "state")
  )

  # SEs must be identical: N=39 with 39 clusters = 1 obs per cluster = HC1
  expect_equal(res_robust$se_att, res_cluster$se_att, tolerance = 1e-10)
  # ATT must also be identical (same point estimate)
  expect_equal(res_robust$att, res_cluster$att, tolerance = 1e-12)
})
