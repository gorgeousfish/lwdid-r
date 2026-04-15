# test-dispatch-estimator.R — dispatch_estimator() unit tests
# Story E6-05, TC-6.5.1 to TC-6.5.46

make_test_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  d <- rbinom(n, 1, plogis(0.5 * x1 - 0.3 * x2))
  y <- 1 + 2 * d + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  data.frame(y = y, d = d, x1 = x1, x2 = x2)
}

# === Group 1: Routing correctness (TC-6.5.1 to TC-6.5.4) ===
test_that("TC-6.5.1: estimator='ra' routes correctly", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ra")
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
  expect_true(res$ci_lower < res$ci_upper)
  expect_identical(res$estimator, "ra")
  expect_identical(res$inference_dist, "t")
})
test_that("TC-6.5.2: estimator='ipw' routes correctly", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipw")
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
  expect_true(res$ci_lower < res$ci_upper)
  expect_identical(res$estimator, "ipw")
  expect_identical(res$inference_dist, "normal")
})
test_that("TC-6.5.3: estimator='ipwra' routes correctly", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipwra")
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
  expect_true(res$ci_lower < res$ci_upper)
  expect_identical(res$estimator, "ipwra")
  expect_identical(res$inference_dist, "normal")
})
test_that("TC-6.5.4: estimator='psm' routes correctly", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm")
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
  expect_true(res$ci_lower < res$ci_upper)
  expect_identical(res$estimator, "psm")
  expect_identical(res$inference_dist, "normal")
})
# === Group 2: Parameter validation (TC-6.5.5 to TC-6.5.18) ===
test_that("TC-6.5.5: ipwra + controls=NULL errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = NULL, estimator = "ipwra"),
    class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.6: ipw + controls=NULL + ps_controls=NULL errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = NULL, ps_controls = NULL, estimator = "ipw"),
    class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.7: psm + controls=NULL + ps_controls=NULL errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = NULL, ps_controls = NULL, estimator = "psm"),
    class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.8: trim_threshold=0 (IPW) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipw",
    trim_threshold = 0), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.9: trim_threshold=0.5 (IPW) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipw",
    trim_threshold = 0.5), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.10: trim_threshold=-0.1 (IPWRA) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipwra",
    trim_threshold = -0.1), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.11: n_neighbors=0 (PSM) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    n_neighbors = 0), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.12: n_neighbors=1.5 (PSM) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    n_neighbors = 1.5), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.13: n_neighbors=-1 (PSM) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    n_neighbors = -1), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.14: caliper=-1 (PSM) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    caliper = -1), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.15: caliper=Inf (PSM) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    caliper = Inf), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.16: caliper=NaN (PSM) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    caliper = NaN), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.17: with_replacement='yes' (PSM) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    with_replacement = "yes"), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.18: se_method='invalid' (IPW) errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipw",
    se_method = "invalid"), class = "lwdid_invalid_parameter")
})
# === Group 3: match.arg validation (TC-6.5.24 to TC-6.5.25) ===
test_that("TC-6.5.24: match_order='invalid' errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    match_order = "invalid"), class = "error")
})
test_that("TC-6.5.25: caliper_scale='invalid' errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    caliper_scale = "invalid"), class = "error")
})
# === Group 4: se_method valid combinations (TC-6.5.26 to TC-6.5.29) ===
test_that("TC-6.5.26: se_method='analytical' for IPW works", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipw",
    se_method = "analytical")
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
})
test_that("TC-6.5.27: se_method='analytical' for IPWRA works", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipwra",
    se_method = "analytical")
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
})
test_that("TC-6.5.28: se_method='abadie_imbens' for PSM works", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    se_method = "abadie_imbens")
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
})
test_that("TC-6.5.29: se_method='bootstrap' for IPW works", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipw",
    se_method = "bootstrap", boot_reps = 10L, seed = 42)
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
})
# === Group 5: Alpha validation (TC-6.5.39 to TC-6.5.40) ===
test_that("TC-6.5.39: alpha=0 errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ra",
    alpha = 0), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.40: alpha=1 errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ra",
    alpha = 1), class = "lwdid_invalid_parameter")
})
# === Group 6: n_neighbors type handling (TC-6.5.41 to TC-6.5.43) ===
test_that("TC-6.5.41: n_neighbors=1 (numeric) works for PSM", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm", n_neighbors = 1)
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
})
test_that("TC-6.5.42: n_neighbors=1.5 errors for PSM", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm",
    n_neighbors = 1.5), class = "lwdid_invalid_parameter")
})
test_that("TC-6.5.43: n_neighbors=2.0 (numeric) works for PSM", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "psm", n_neighbors = 2.0)
  expect_true(is.numeric(res$att))
  expect_true(res$se > 0)
})
# === Group 7: Controls empty normalization (TC-6.5.44) ===
test_that("TC-6.5.44: controls=character(0) + ipwra errors", {
  df <- make_test_data(n = 20)
  expect_error(dispatch_estimator(df, "y", "d",
    controls = character(0), estimator = "ipwra"),
    class = "lwdid_invalid_parameter")
})
# === Group 8: Return structure (TC-6.5.45 to TC-6.5.46) ===
test_that("TC-6.5.45: all estimators return estimator field", {
  df <- make_test_data()
  for (est in c("ra", "ipw", "ipwra", "psm")) {
    res <- dispatch_estimator(df, "y", "d",
      controls = c("x1", "x2"), estimator = est)
    expect_identical(res$estimator, est,
      info = paste("estimator field for", est))
  }
})
test_that("TC-6.5.46: inference_dist t for RA, normal for others", {
  df <- make_test_data()
  expected <- c(ra = "t", ipw = "normal",
    ipwra = "normal", psm = "normal")
  for (est in names(expected)) {
    res <- dispatch_estimator(df, "y", "d",
      controls = c("x1", "x2"), estimator = est)
    expect_identical(res$inference_dist, expected[[est]],
      info = paste("inference_dist for", est))
  }
})
# === Group 9: Inference distribution verification (TC-6.5.21-23) ===
test_that("TC-6.5.21: IPW CI matches normal distribution", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ipw")
  z <- qnorm(1 - 0.05 / 2)
  expected_lower <- res$att - z * res$se
  expected_upper <- res$att + z * res$se
  expect_equal(res$ci_lower, expected_lower, tolerance = 1e-6)
  expect_equal(res$ci_upper, expected_upper, tolerance = 1e-6)
})
test_that("TC-6.5.22: RA CI wider than normal (uses t-dist)", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ra")
  z <- qnorm(1 - 0.05 / 2)
  ci_width_normal <- 2 * z * res$se
  ci_width_actual <- res$ci_upper - res$ci_lower
  expect_true(ci_width_actual > ci_width_normal,
    info = "t-dist CI should be wider than normal CI")
})
test_that("TC-6.5.23: RA CI consistent with t-distribution", {
  df <- make_test_data()
  res <- dispatch_estimator(df, "y", "d",
    controls = c("x1", "x2"), estimator = "ra")
  n_treat <- sum(df$d == 1)
  n_ctrl <- sum(df$d == 0)
  n_controls <- 2
  df_t <- n_treat + n_ctrl - n_controls - 2
  t_crit <- qt(1 - 0.05 / 2, df = df_t)
  expected_lower <- res$att - t_crit * res$se
  expected_upper <- res$att + t_crit * res$se
  expect_equal(res$ci_lower, expected_lower, tolerance = 1e-4)
  expect_equal(res$ci_upper, expected_upper, tolerance = 1e-4)
})
