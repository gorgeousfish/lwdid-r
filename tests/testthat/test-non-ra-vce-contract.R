# ============================================================================
# test-non-ra-vce-contract.R
# Non-RA estimators must not accept RA/OLS VCE options.
# ============================================================================

make_non_ra_vce_ct_data <- function(n_units = 40L, n_periods = 8L) {
  df <- expand.grid(id = seq_len(n_units), year = 2000L + seq_len(n_periods))
  df <- df[order(df$id, df$year), ]

  set.seed(20260525)
  unit_x <- stats::rnorm(n_units)
  df$d <- as.integer(df$id <= n_units / 2L)
  df$post <- as.integer(df$year >= 2005L)
  df$x <- unit_x[df$id]
  df$cluster <- rep(seq_len(10L), each = 4L)[df$id]
  df$y <- 0.2 * df$x + 0.5 * df$d * df$post + stats::rnorm(nrow(df))
  rownames(df) <- NULL
  df
}

test_that("Castle staggered non-RA rejects cluster vce before row labeling", {
  data("castle", package = "lwdid")

  expect_error(
    lwdid(
      data = castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      estimator = "ipw",
      aggregate = "none",
      control_group = "not_yet_treated",
      controls = c("income", "unemployrt", "poverty"),
      vce = "cluster",
      cluster_var = "sid",
      verbose = "quiet"
    ),
    regexp = "vce is only supported when estimator='ra'",
    class = "lwdid_invalid_parameter"
  )
})

test_that("common-timing non-RA rejects bootstrap vce before WCB routing", {
  df <- make_non_ra_vce_ct_data()

  expect_error(
    lwdid(
      data = df,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      estimator = "ipw",
      controls = "x",
      vce = "bootstrap",
      cluster_var = "cluster",
      verbose = "quiet"
    ),
    regexp = "vce is only supported when estimator='ra'",
    class = "lwdid_invalid_parameter"
  )
})

test_that("common-timing non-RA still accepts estimator-specific bootstrap SE", {
  df <- make_non_ra_vce_ct_data(n_units = 80L)

  fit <- suppressWarnings(lwdid(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    estimator = "ipw",
    controls = "x",
    se_method = "bootstrap",
    boot_reps = 10L,
    seed = 123L,
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  expect_s3_class(fit, "lwdid_result")
  expect_identical(fit$estimator, "ipw")
  expect_identical(fit$vce_type, "bootstrap")
  expect_true(is.finite(fit$se_att))
  expect_false(is.null(fit$propensity_scores))
  expect_length(fit$propensity_scores, 80L)
  expect_true(is.finite(fit$weights_cv))
  propensity_diag <- get_diagnostics(fit, "propensity")
  expect_type(propensity_diag, "list")
  expect_identical(propensity_diag$estimator, "ipw")
  expect_identical(propensity_diag$scope, "common_timing")
  expect_true(is.finite(propensity_diag$ps_min))
  expect_true(is.finite(propensity_diag$ps_max))
  expect_match(summary(fit)$vce_description, "Bootstrap")
  expect_match(summary(fit)$diagnostics_summary$propensity, "ps_range")
})
