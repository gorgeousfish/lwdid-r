# ============================================================================
# Tests for compute_inference()
# Story E3-04: Unified t-Inference for ATT
#
# Test Groups:
#   Group 1: Basic t-inference correctness (E3-04.2)
#
# compute_inference(att, vcov_mat, df, alpha = 0.05, coef_index = 2L)
# Returns: list(att, se, t_stat, df, pvalue, ci_lower, ci_upper)
#
# Formulas:
#   se       = sqrt(vcov_mat[coef_index, coef_index])
#   t_stat   = att / se
#   pvalue   = 2 * pt(-abs(t_stat), df)
#   ci_lower = att - qt(1 - alpha/2, df) * se
#   ci_upper = att + qt(1 - alpha/2, df) * se
# ============================================================================

# ============================================================================
# Group 1: Basic t-inference correctness (E3-04.2)
# ============================================================================

test_that("known VCE matrix produces correct se, t_stat, pvalue, CI", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  att <- 2.5
  df_val <- 10L
  coef_index <- 2L

  res <- compute_inference(att, vcov_mat, df = df_val, coef_index = coef_index)

  # Expected values computed from R's own pt/qt
  expected_se <- sqrt(0.25)  # 0.5
  expected_t  <- att / expected_se  # 5.0
  expected_p  <- 2 * pt(-abs(expected_t), df_val)
  expected_ci_lower <- att - qt(0.975, df_val) * expected_se
  expected_ci_upper <- att + qt(0.975, df_val) * expected_se

  expect_equal(res$se, expected_se, tolerance = 1e-12)
  expect_equal(res$t_stat, expected_t, tolerance = 1e-12)
  expect_equal(res$pvalue, expected_p, tolerance = 1e-12)
  expect_equal(res$ci_lower, expected_ci_lower, tolerance = 1e-12)
  expect_equal(res$ci_upper, expected_ci_upper, tolerance = 1e-12)
})

test_that("return structure is list with exactly 7 named numeric scalars", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  res <- compute_inference(att = 1.0, vcov_mat = vcov_mat, df = 10)

  expect_true(is.list(res))

  expected_names <- c("att", "se", "t_stat", "df", "pvalue",
                       "ci_lower", "ci_upper")
  expect_identical(sort(names(res)), sort(expected_names))
  expect_equal(length(res), 7L)

  # Each element is a numeric scalar (length-1 numeric vector)
  for (nm in expected_names) {
    expect_true(is.numeric(res[[nm]]),
                info = paste(nm, "should be numeric"))
    expect_equal(length(res[[nm]]), 1L,
                 info = paste(nm, "should be scalar"))
  }
})

test_that("att = 0 yields t_stat = 0, pvalue = 1, symmetric CI about 0", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  res <- compute_inference(att = 0, vcov_mat = vcov_mat, df = 10)

  expect_equal(res$t_stat, 0, tolerance = 1e-12)
  expect_equal(res$pvalue, 1.0, tolerance = 1e-12)
  # CI symmetric about 0: ci_lower = -ci_upper
  expect_equal(res$ci_lower, -res$ci_upper, tolerance = 1e-12)
})

test_that("negative ATT produces negative t_stat, same pvalue as positive, CI below zero", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  df_val <- 10L

  res_neg <- compute_inference(att = -3.0, vcov_mat = vcov_mat, df = df_val)
  res_pos <- compute_inference(att =  3.0, vcov_mat = vcov_mat, df = df_val)

  # t_stat is negative

  expect_true(res_neg$t_stat < 0)

  # pvalue same as for positive att of same magnitude
  expect_equal(res_neg$pvalue, res_pos$pvalue, tolerance = 1e-12)

  # CI entirely below zero
  expect_true(res_neg$ci_upper < 0)
})

test_that("large t_stat yields pvalue essentially zero", {
  # att = 100, se = sqrt(0.01) = 0.1 -> t_stat = 1000
  vcov_mat <- matrix(c(1, 0, 0, 0.01), nrow = 2)
  res <- compute_inference(att = 100, vcov_mat = vcov_mat, df = 50)

  expect_equal(res$t_stat, 1000, tolerance = 1e-12)
  expect_true(res$pvalue < 1e-100)
})

test_that("coef_index = 1 extracts variance from first diagonal element", {
  # First diagonal = 0.16 (se = 0.4), second diagonal = 0.25 (se = 0.5)
  vcov_mat <- matrix(c(0.16, 0, 0, 0.25), nrow = 2)
  att <- 2.0
  df_val <- 20L

  res_idx1 <- compute_inference(att, vcov_mat, df = df_val, coef_index = 1L)
  res_idx2 <- compute_inference(att, vcov_mat, df = df_val, coef_index = 2L)

  # coef_index = 1 should use var = 0.16 -> se = 0.4
  expect_equal(res_idx1$se, sqrt(0.16), tolerance = 1e-12)
  # coef_index = 2 should use var = 0.25 -> se = 0.5
  expect_equal(res_idx2$se, sqrt(0.25), tolerance = 1e-12)

  # Different SEs -> different t_stats
  expect_equal(res_idx1$t_stat, att / sqrt(0.16), tolerance = 1e-12)
  expect_equal(res_idx2$t_stat, att / sqrt(0.25), tolerance = 1e-12)
})

test_that("different df values produce different CI widths and pvalues", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  att <- 2.0

  res_df5   <- compute_inference(att, vcov_mat, df = 5)
  res_df100 <- compute_inference(att, vcov_mat, df = 100)

  # CI should be wider for df = 5 (heavier tails)
  ci_width_5   <- res_df5$ci_upper   - res_df5$ci_lower
  ci_width_100 <- res_df100$ci_upper - res_df100$ci_lower
  expect_true(ci_width_5 > ci_width_100)

  # pvalue should be larger for df = 5 (heavier tails)
  expect_true(res_df5$pvalue > res_df100$pvalue)
})

test_that("att and df are passed through exactly", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  att_in <- 3.14159
  df_in  <- 42L

  res <- compute_inference(att = att_in, vcov_mat = vcov_mat, df = df_in)

  expect_identical(res$att, att_in)
  expect_identical(res$df, df_in)
})

# ============================================================================
# Group 2: t-distribution vs normal distribution distinction (Task E3-04.3)
# ============================================================================

test_that("compute_inference: t-distribution vs normal", {

  # --------------------------------------------------------------------------
  # Test 1: Small df (df=3) — t vs normal p-value difference
  # --------------------------------------------------------------------------
  vcov_1 <- matrix(c(1, 0, 0, 1.0), nrow = 2)
  res_1 <- compute_inference(att = 2.0, vcov_mat = vcov_1, df = 3,
                             coef_index = 2L)

  t_stat_1 <- 2.0  # att / se = 2.0 / 1.0
  t_pvalue_1 <- 2 * pt(-2.0, 3)
  norm_pvalue_1 <- 2 * pnorm(-2.0)

  # compute_inference() must use t-distribution

  expect_equal(res_1$pvalue, t_pvalue_1, tolerance = 1e-12)

  # t-distribution has heavier tails => larger p-value than normal
  expect_true(t_pvalue_1 > norm_pvalue_1)

  # The difference should be substantial for df=3
  diff_1 <- t_pvalue_1 - norm_pvalue_1
  expect_true(diff_1 > 0.01)

  # --------------------------------------------------------------------------
  # Test 2: df=5 — another small df test
  # --------------------------------------------------------------------------
  vcov_2 <- matrix(c(1, 0, 0, 1.0), nrow = 2)
  res_2 <- compute_inference(att = 3.0, vcov_mat = vcov_2, df = 5,
                             coef_index = 2L)

  t_pvalue_2 <- 2 * pt(-3.0, 5)
  norm_pvalue_2 <- 2 * pnorm(-3.0)

  # Must match t-distribution exactly
  expect_equal(res_2$pvalue, t_pvalue_2, tolerance = 1e-12)

  # t p-value > normal p-value
  expect_true(t_pvalue_2 > norm_pvalue_2)

  # --------------------------------------------------------------------------
  # Test 3: df=1 (Cauchy distribution) — extreme case
  # --------------------------------------------------------------------------
  vcov_3 <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  res_3 <- compute_inference(att = 1.0, vcov_mat = vcov_3, df = 1,
                             coef_index = 2L)

  # Should not crash — normal computation
  expect_true(is.numeric(res_3$pvalue))
  expect_true(is.numeric(res_3$ci_lower))
  expect_true(is.numeric(res_3$ci_upper))

  # qt(0.975, 1) is the famous Cauchy critical value ≈ 12.706
  expect_equal(qt(0.975, 1), 12.706, tolerance = 1e-3)

  # CI should be extremely wide: 2 * 12.706 * 0.5 ≈ 12.706
  ci_width_3 <- res_3$ci_upper - res_3$ci_lower
  expect_true(ci_width_3 > 12)

  # t p-value >> normal p-value for df=1
  t_pvalue_3 <- res_3$pvalue
  # t_stat = 1.0 / 0.5 = 2.0
  norm_pvalue_3 <- 2 * pnorm(-2.0)
  expect_true(t_pvalue_3 > norm_pvalue_3)

  # --------------------------------------------------------------------------
  # Test 4: Large df (df=1000) — t approaches normal but not equal
  # --------------------------------------------------------------------------
  vcov_4 <- matrix(c(1, 0, 0, 1.0), nrow = 2)
  res_4 <- compute_inference(att = 2.0, vcov_mat = vcov_4, df = 1000,
                             coef_index = 2L)

  t_pvalue_4 <- 2 * pt(-2.0, 1000)
  norm_pvalue_4 <- 2 * pnorm(-2.0)

  # They should be close (difference < 1e-3)
  expect_true(abs(t_pvalue_4 - norm_pvalue_4) < 1e-3)

  # But NOT exactly equal
  expect_true(t_pvalue_4 != norm_pvalue_4)

  # compute_inference() must use t-distribution (matches pt, not pnorm)
  expect_equal(res_4$pvalue, t_pvalue_4, tolerance = 1e-12)

  # --------------------------------------------------------------------------
  # Test 5: df=10000 — near-asymptotic
  # --------------------------------------------------------------------------
  vcov_5 <- matrix(c(1, 0, 0, 1.0), nrow = 2)
  res_5 <- compute_inference(att = 1.96, vcov_mat = vcov_5, df = 10000,
                             coef_index = 2L)

  t_pvalue_5 <- 2 * pt(-1.96, 10000)
  norm_pvalue_5 <- 2 * pnorm(-1.96)

  # Very close to normal (difference < 1e-4)
  expect_true(abs(t_pvalue_5 - norm_pvalue_5) < 1e-4)

  # But still verify compute_inference() matches pt() exactly
  expect_equal(res_5$pvalue, t_pvalue_5, tolerance = 1e-12)

  # --------------------------------------------------------------------------
  # Test 6: CI width comparison across df values
  # --------------------------------------------------------------------------
  alpha <- 0.05
  vcov_ci <- matrix(c(1, 0, 0, 1.0), nrow = 2)  # se = 1.0

  res_df3 <- compute_inference(att = 1.0, vcov_mat = vcov_ci, df = 3,
                               alpha = alpha, coef_index = 2L)
  res_df10 <- compute_inference(att = 1.0, vcov_mat = vcov_ci, df = 10,
                                alpha = alpha, coef_index = 2L)
  res_df100 <- compute_inference(att = 1.0, vcov_mat = vcov_ci, df = 100,
                                 alpha = alpha, coef_index = 2L)
  res_df10k <- compute_inference(att = 1.0, vcov_mat = vcov_ci, df = 10000,
                                 alpha = alpha, coef_index = 2L)

  width_df3 <- res_df3$ci_upper - res_df3$ci_lower
  width_df10 <- res_df10$ci_upper - res_df10$ci_lower
  width_df100 <- res_df100$ci_upper - res_df100$ci_lower
  width_df10k <- res_df10k$ci_upper - res_df10k$ci_lower

  # CI width must be monotonically decreasing as df increases
  expect_true(width_df3 > width_df10)
  expect_true(width_df10 > width_df100)
  expect_true(width_df100 > width_df10k)

  # All t-based CI widths must exceed the normal CI width
  normal_ci_width <- 2 * qnorm(0.975) * 1.0
  expect_true(width_df3 > normal_ci_width)
  expect_true(width_df10 > normal_ci_width)
  expect_true(width_df100 > normal_ci_width)
  expect_true(width_df10k > normal_ci_width)
})

# ============================================================================
# Group 3: Alpha parameter tests (Task E3-04.4)
# ============================================================================

test_that("compute_inference: alpha controls CI width correctly", {
  att <- 2.0
  vcov_mat <- matrix(c(1, 0, 0, 1.0), nrow = 2)  # se = 1.0
  df_val <- 20L
  se <- 1.0


  # --- Test 1: alpha=0.05 (95% CI) ---
  res_05 <- compute_inference(att, vcov_mat, df = df_val, alpha = 0.05,
                              coef_index = 2L)
  t_crit_05 <- qt(0.975, df_val)
  expect_equal(res_05$ci_lower, att - t_crit_05 * se, tolerance = 1e-12)
  expect_equal(res_05$ci_upper, att + t_crit_05 * se, tolerance = 1e-12)

  # --- Test 2: alpha=0.10 (90% CI) ---
  res_10 <- compute_inference(att, vcov_mat, df = df_val, alpha = 0.10,
                              coef_index = 2L)
  t_crit_10 <- qt(0.95, df_val)
  expect_equal(res_10$ci_lower, att - t_crit_10 * se, tolerance = 1e-12)
  expect_equal(res_10$ci_upper, att + t_crit_10 * se, tolerance = 1e-12)

  # 90% CI must be NARROWER than 95% CI
  width_90 <- res_10$ci_upper - res_10$ci_lower
  width_95 <- res_05$ci_upper - res_05$ci_lower
  expect_true(width_90 < width_95)

  # --- Test 3: alpha=0.01 (99% CI) ---
  res_01 <- compute_inference(att, vcov_mat, df = df_val, alpha = 0.01,
                              coef_index = 2L)
  t_crit_01 <- qt(0.995, df_val)
  expect_equal(res_01$ci_lower, att - t_crit_01 * se, tolerance = 1e-12)
  expect_equal(res_01$ci_upper, att + t_crit_01 * se, tolerance = 1e-12)

  # 99% CI must be WIDER than 95% CI
  width_99 <- res_01$ci_upper - res_01$ci_lower
  expect_true(width_99 > width_95)

  # --- Test 4: Width ordering: 99% > 95% > 90% (strict) ---
  expect_true(width_99 > width_95)
  expect_true(width_95 > width_90)

  # --- Test 5: p-value and t_stat unchanged by alpha ---
  expect_equal(res_05$pvalue, res_10$pvalue, tolerance = 1e-12)
  expect_equal(res_05$pvalue, res_01$pvalue, tolerance = 1e-12)
  expect_equal(res_05$t_stat, res_10$t_stat, tolerance = 1e-12)
  expect_equal(res_05$t_stat, res_01$t_stat, tolerance = 1e-12)
})

test_that("compute_inference: invalid alpha values throw lwdid_invalid_parameter", {
  vcov_mat <- matrix(c(1, 0, 0, 1), nrow = 2)

  # alpha = 0 (boundary, not in open interval)
  expect_error(
    compute_inference(att = 1, vcov_mat = vcov_mat, df = 10, alpha = 0),
    class = "lwdid_invalid_parameter"
  )

  # alpha = 1 (boundary, not in open interval)
  expect_error(
    compute_inference(att = 1, vcov_mat = vcov_mat, df = 10, alpha = 1),
    class = "lwdid_invalid_parameter"
  )

  # alpha = NA
  expect_error(
    compute_inference(att = 1, vcov_mat = vcov_mat, df = 10, alpha = NA),
    class = "lwdid_invalid_parameter"
  )

  # alpha = "0.05" (string, not numeric)
  expect_error(
    compute_inference(att = 1, vcov_mat = vcov_mat, df = 10, alpha = "0.05"),
    class = "lwdid_invalid_parameter"
  )

  # alpha = c(0.05, 0.10) (vector length > 1)
  expect_error(
    compute_inference(att = 1, vcov_mat = vcov_mat, df = 10,
                      alpha = c(0.05, 0.10)),
    class = "lwdid_invalid_parameter"
  )

  # alpha = -0.5 (negative)
  expect_error(
    compute_inference(att = 1, vcov_mat = vcov_mat, df = 10, alpha = -0.5),
    class = "lwdid_invalid_parameter"
  )

  # alpha = 1.5 (> 1)
  expect_error(
    compute_inference(att = 1, vcov_mat = vcov_mat, df = 10, alpha = 1.5),
    class = "lwdid_invalid_parameter"
  )

  # alpha = NaN
  expect_error(
    compute_inference(att = 1, vcov_mat = vcov_mat, df = 10, alpha = NaN),
    class = "lwdid_invalid_parameter"
  )
})

# ============================================================================
# Group 3: Alpha parameter tests (E3-04.4)
# ============================================================================

test_that("alpha=0.05 (95% CI) uses qt(0.975, df)", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  att <- 2.0
  df_val <- 10L
  se <- sqrt(0.25)  # 0.5

  res <- compute_inference(att, vcov_mat, df = df_val, alpha = 0.05)

  t_crit <- qt(0.975, df_val)
  expect_equal(res$ci_lower, att - t_crit * se, tolerance = 1e-12)
  expect_equal(res$ci_upper, att + t_crit * se, tolerance = 1e-12)
})

test_that("alpha=0.10 (90% CI) uses qt(0.95, df) and is narrower than 95%", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  att <- 2.0
  df_val <- 10L
  se <- sqrt(0.25)

  res_90 <- compute_inference(att, vcov_mat, df = df_val, alpha = 0.10)
  res_95 <- compute_inference(att, vcov_mat, df = df_val, alpha = 0.05)

  t_crit_90 <- qt(0.95, df_val)
  expect_equal(res_90$ci_lower, att - t_crit_90 * se, tolerance = 1e-12)
  expect_equal(res_90$ci_upper, att + t_crit_90 * se, tolerance = 1e-12)

  # 90% CI must be narrower than 95% CI
  width_90 <- res_90$ci_upper - res_90$ci_lower
  width_95 <- res_95$ci_upper - res_95$ci_lower
  expect_true(width_90 < width_95)
})

test_that("alpha=0.01 (99% CI) uses qt(0.995, df) and is wider than 95%", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  att <- 2.0
  df_val <- 10L
  se <- sqrt(0.25)

  res_99 <- compute_inference(att, vcov_mat, df = df_val, alpha = 0.01)
  res_95 <- compute_inference(att, vcov_mat, df = df_val, alpha = 0.05)

  t_crit_99 <- qt(0.995, df_val)
  expect_equal(res_99$ci_lower, att - t_crit_99 * se, tolerance = 1e-12)
  expect_equal(res_99$ci_upper, att + t_crit_99 * se, tolerance = 1e-12)

  # 99% CI must be wider than 95% CI
  width_99 <- res_99$ci_upper - res_99$ci_lower
  width_95 <- res_95$ci_upper - res_95$ci_lower
  expect_true(width_99 > width_95)
})

test_that("alpha=0 throws lwdid_invalid_parameter", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, 10, alpha = 0),
    class = "lwdid_invalid_parameter"
  )
})

test_that("alpha=1 throws lwdid_invalid_parameter", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, 10, alpha = 1),
    class = "lwdid_invalid_parameter"
  )
})

test_that("alpha=NA throws lwdid_invalid_parameter", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, 10, alpha = NA),
    class = "lwdid_invalid_parameter"
  )
})

test_that("alpha='0.05' (non-numeric) throws lwdid_invalid_parameter", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, 10, alpha = "0.05"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("alpha=c(0.05, 0.10) (length > 1) throws lwdid_invalid_parameter", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, 10, alpha = c(0.05, 0.10)),
    class = "lwdid_invalid_parameter"
  )
})

# ============================================================================
# Group 4: Input validation tests (E3-04.5)
# ============================================================================

test_that("att = NA_real_ throws lwdid_numerical", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(NA_real_, vcov_mat, 10),
    class = "lwdid_numerical"
  )
})

test_that("att = NaN throws lwdid_numerical", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(NaN, vcov_mat, 10),
    class = "lwdid_numerical"
  )
})

test_that("coef_index = 0L throws lwdid_invalid_parameter", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, 10, coef_index = 0L),
    class = "lwdid_invalid_parameter"
  )
})

test_that("coef_index = 3L (out of bounds for 2x2) throws lwdid_invalid_parameter", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, 10, coef_index = 3L),
    class = "lwdid_invalid_parameter"
  )
})

test_that("df = 0 throws lwdid_insufficient_data", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, df = 0),
    class = "lwdid_insufficient_data"
  )
})

test_that("df = -1 throws lwdid_insufficient_data", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, df = -1),
    class = "lwdid_insufficient_data"
  )
})

test_that("df = NA throws lwdid_insufficient_data", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  expect_error(
    compute_inference(1, vcov_mat, df = NA),
    class = "lwdid_insufficient_data"
  )
})

test_that("VCE diagonal NA throws lwdid_numerical", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  vcov_mat[2, 2] <- NA
  expect_error(
    compute_inference(1, vcov_mat, 10),
    class = "lwdid_numerical"
  )
})

test_that("VCE diagonal NaN throws lwdid_numerical", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  vcov_mat[2, 2] <- NaN
  expect_error(
    compute_inference(1, vcov_mat, 10),
    class = "lwdid_numerical"
  )
})

test_that("VCE diagonal negative (-1e-15) throws lwdid_numerical", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  vcov_mat[2, 2] <- -1e-15
  expect_error(
    compute_inference(1, vcov_mat, 10),
    class = "lwdid_numerical"
  )
})

test_that("validation order: alpha checked before att", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  # Both alpha=2 (invalid) and att=NA (invalid) — alpha error fires first
  expect_error(
    compute_inference(att = NA, vcov_mat, 10, alpha = 2),
    class = "lwdid_invalid_parameter"
  )
})

test_that("validation order: att checked before df", {
  vcov_mat <- matrix(c(1, 0, 0, 1), nrow = 2)
  # Both att=NA (invalid) and df=-1 (invalid) — att error fires first
  expect_error(
    compute_inference(att = NA, vcov_mat = vcov_mat, df = -1, alpha = 0.05),
    class = "lwdid_numerical"
  )
})

test_that("validation order: coef_index checked before df", {
  vcov_mat <- matrix(c(1, 0, 0, 1), nrow = 2)
  # Both coef_index=5L (invalid) and df=-1 (invalid) — coef_index fires first
  expect_error(
    compute_inference(
      att = 1.0, vcov_mat = vcov_mat, df = -1,
      alpha = 0.05, coef_index = 5L
    ),
    class = "lwdid_invalid_parameter"
  )
})


# ============================================================================
# Group 5: Boundary case tests (E3-04.6)
# Design.md §6 boundary cases BC-01 through BC-27
# ============================================================================

test_that("BC-15: SE=0 with att!=0 — Inf t_stat, pvalue=0, CI=[att,att]", {
  vcov_mat <- matrix(c(1, 0, 0, 0), nrow = 2)
  att <- 5.0
  df_val <- 10L

  expect_warning(
    res <- compute_inference(att, vcov_mat, df = df_val, coef_index = 2L),
    class = "lwdid_numerical"
  )

  expect_equal(res$t_stat, Inf)
  expect_equal(res$pvalue, 0)
  expect_equal(res$ci_lower, 5.0)
  expect_equal(res$ci_upper, 5.0)
  expect_equal(res$se, 0)
})

test_that("BC-16: SE=0 with att=0 — NaN t_stat, NaN pvalue, CI=[0,0]", {
  vcov_mat <- matrix(c(1, 0, 0, 0), nrow = 2)
  att <- 0
  df_val <- 10L

  expect_warning(
    res <- compute_inference(att, vcov_mat, df = df_val, coef_index = 2L),
    class = "lwdid_numerical"
  )

  expect_true(is.nan(res$t_stat))
  expect_true(is.nan(res$pvalue))
  expect_equal(res$ci_lower, 0)
  expect_equal(res$ci_upper, 0)
})

test_that("BC-11: att=0 with se>0 — t_stat=0, pvalue=1, symmetric CI", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  att <- 0
  df_val <- 10L

  res <- compute_inference(att, vcov_mat, df = df_val, coef_index = 2L)

  expect_equal(res$t_stat, 0, tolerance = 1e-12)
  expect_equal(res$pvalue, 1.0, tolerance = 1e-12)
  expect_equal(res$ci_lower, -res$ci_upper, tolerance = 1e-12)
})

test_that("BC-01: df=1 Cauchy — extremely wide CI", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  att <- 1.0
  df_val <- 1L

  res <- compute_inference(att, vcov_mat, df = df_val, coef_index = 2L)

  # qt(0.975, 1) is the Cauchy critical value
  expect_equal(qt(0.975, 1), 12.706, tolerance = 1e-3)

  # CI width = 2 * qt(0.975, 1) * se = 2 * 12.706 * 0.5
  ci_width <- res$ci_upper - res$ci_lower
  expect_true(ci_width > 12)

  # All results are finite numeric
  expect_true(is.finite(res$t_stat))
  expect_true(is.finite(res$pvalue))
  expect_true(is.finite(res$ci_lower))
  expect_true(is.finite(res$ci_upper))
})

test_that("BC-05: large df=10000 — t approaches normal but not equal", {
  vcov_mat <- matrix(c(1, 0, 0, 1.0), nrow = 2)
  att <- 2.0
  df_val <- 10000L

  res <- compute_inference(att, vcov_mat, df = df_val, coef_index = 2L)

  # p-value close to normal p-value (within 1e-4)
  norm_pvalue <- 2 * pnorm(-2.0)
  expect_true(abs(res$pvalue - norm_pvalue) < 1e-4)

  # But matches pt() exactly (tolerance 1e-12)
  t_pvalue <- 2 * pt(-abs(att / 1.0), df_val)
  expect_equal(res$pvalue, t_pvalue, tolerance = 1e-12)

  # CI width close to normal CI width but not identical
  ci_width <- res$ci_upper - res$ci_lower
  normal_ci_width <- 2 * qnorm(0.975) * 1.0
  expect_true(abs(ci_width - normal_ci_width) < 0.01)
  expect_false(identical(ci_width, normal_ci_width))
})

test_that("BC-14: extreme att=1e10 — huge t_stat, pvalue near 0", {
  vcov_mat <- matrix(c(1, 0, 0, 1.0), nrow = 2)
  att <- 1e10
  df_val <- 50L

  res <- compute_inference(att, vcov_mat, df = df_val, coef_index = 2L)

  expect_equal(res$t_stat, 1e10, tolerance = 1e-12)
  expect_true(res$pvalue < 1e-100)
  expect_true(res$ci_lower > 0)
})

test_that("BC-17: tiny SE, var_att=1e-25 — finite results", {
  vcov_mat <- matrix(c(1, 0, 0, 1e-25), nrow = 2)
  att <- 1e-10
  df_val <- 30L

  res <- compute_inference(att, vcov_mat, df = df_val, coef_index = 2L)

  # se = sqrt(1e-25)
  expected_se <- sqrt(1e-25)
  expect_equal(res$se, expected_se, tolerance = 1e-20)

  # t_stat = att / se is finite
  expect_true(is.finite(res$t_stat))
  expected_t <- att / expected_se
  expect_equal(res$t_stat, expected_t, tolerance = 1e-6)

  # pvalue is finite and in [0, 1]
  expect_true(is.finite(res$pvalue))
  expect_true(res$pvalue >= 0 && res$pvalue <= 1)

  # CI bounds are finite
  expect_true(is.finite(res$ci_lower))
  expect_true(is.finite(res$ci_upper))
})

# ============================================================================
# Group 6: R-Python numerical consistency (E3-04.7)
# ============================================================================

test_that("p-value matches 2*pt(-abs(t_stat), df) for multiple combos", {
  # Helper: build vcov so that se = 1.0, then att = desired t_stat
  make_vcov <- function() matrix(c(1, 0, 0, 1.0), nrow = 2)

  combos <- list(
    list(t_stat = 2.5,  df = 10L),
    list(t_stat = 1.96, df = 30L),
    list(t_stat = 3.0,  df = 5L),
    list(t_stat = 0.5,  df = 100L)
  )

  for (combo in combos) {
    # With se = 1.0, att = t_stat gives the desired t_stat
    att <- combo$t_stat
    df_val <- combo$df
    expected_p <- 2 * pt(-abs(combo$t_stat), df_val)

    res <- compute_inference(att, make_vcov(), df = df_val)

    expect_equal(
      res$pvalue, expected_p,
      tolerance = 1e-12,
      info = sprintf("t_stat=%.2f, df=%d", combo$t_stat, df_val)
    )
  }
})

test_that("CI bounds match att +/- qt(0.975, df) * se for multiple df", {
  att <- 2.0
  se <- 0.5  # var = 0.25
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)

  df_vals <- c(10L, 30L, 120L)

  for (df_val in df_vals) {
    t_crit <- qt(0.975, df_val)
    expected_lower <- att - t_crit * se
    expected_upper <- att + t_crit * se

    res <- compute_inference(att, vcov_mat, df = df_val)

    expect_equal(
      res$ci_lower, expected_lower,
      tolerance = 1e-12,
      info = sprintf("df=%d ci_lower", df_val)
    )
    expect_equal(
      res$ci_upper, expected_upper,
      tolerance = 1e-12,
      info = sprintf("df=%d ci_upper", df_val)
    )
  }
})

test_that("cross-validation with hardcoded reference values", {
  # Reference values computed from R's pt()/qt() and cross-checked
  # against Python scipy.stats.t (both use incomplete beta function,
  # agree to machine precision ~1e-15).
  #
  # p-value references: 2 * pt(-|t|, df)
  expect_equal(2 * pt(-2.5, 10),  0.031446844236609, tolerance = 1e-12)
  expect_equal(2 * pt(-1.96, 30), 0.059342312896051, tolerance = 1e-12)
  expect_equal(2 * pt(-3.0, 5),   0.030099247897463, tolerance = 1e-12)
  expect_equal(2 * pt(-0.5, 100), 0.618173565830887, tolerance = 1e-12)

  # Critical value references: qt(0.975, df)
  expect_equal(qt(0.975, 10),  2.228138851986274, tolerance = 1e-12)
  expect_equal(qt(0.975, 30),  2.042272456301238, tolerance = 1e-12)
  expect_equal(qt(0.975, 120), 1.979930405082422, tolerance = 1e-12)

  # Verify compute_inference() produces results matching these references
  vcov_mat <- matrix(c(1, 0, 0, 1.0), nrow = 2)  # se = 1.0

  # t_stat=2.5, df=10
  res1 <- compute_inference(att = 2.5, vcov_mat, df = 10L)
  expect_equal(res1$pvalue, 0.031446844236609, tolerance = 1e-12)

  # t_stat=1.96, df=30
  res2 <- compute_inference(att = 1.96, vcov_mat, df = 30L)
  expect_equal(res2$pvalue, 0.059342312896051, tolerance = 1e-12)

  # t_stat=3.0, df=5
  res3 <- compute_inference(att = 3.0, vcov_mat, df = 5L)
  expect_equal(res3$pvalue, 0.030099247897463, tolerance = 1e-12)

  # t_stat=0.5, df=100
  res4 <- compute_inference(att = 0.5, vcov_mat, df = 100L)
  expect_equal(res4$pvalue, 0.618173565830887, tolerance = 1e-12)
})
