# ============================================================================
# test-vm-pretreatment-joint.R
# vibe-math MCP numerical cross-validation for test_parallel_trends_joint()
# Story E7-05: Parallel Trends Joint Test
# Paper: LW2025 Equation (2.3) Parallel Trends Hypothesis
#
# Each verification item uses independently computed expected values
# from vibe-math MCP and R reference calculations.
# ============================================================================

# Helper: build pre_treatment_effects data.frame
make_vm_effects <- function(t_stats, n_treated, n_control,
                             event_time = NULL) {
  K <- length(t_stats)
  att <- t_stats  # se = 1 → att = t_stat
  se <- rep(1.0, K)
  pvalue <- 2 * stats::pt(-abs(t_stats), df = 48)
  if (is.null(event_time)) event_time <- seq(-K, -1L)
  data.frame(
    att = att, se = se, t_stat = t_stats, pvalue = pvalue,
    is_anchor = rep(FALSE, K),
    n_treated = as.integer(n_treated),
    n_control = as.integer(n_control),
    event_time = as.integer(event_time),
    cohort = rep(NA_integer_, K),
    period = seq_len(K),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# VM-7.5.1: F statistic verification
# Input: t_stats = [2.0, 1.5, 0.5], K = 3
# Expected: F = (4 + 2.25 + 0.25) / 3 = 13/6 ≈ 2.16667
# vibe-math result: 2.1666666666666665
# ============================================================================
test_that("VM-7.5.1: F statistic = sum(t²)/K verified by vibe-math", {
  pe <- make_vm_effects(
    t_stats = c(2.0, 1.5, 0.5),
    n_treated = rep(25L, 3), n_control = rep(25L, 3)
  )
  res <- test_parallel_trends_joint(pe, test_type = "f")
  # vibe-math: (2.0^2 + 1.5^2 + 0.5^2) / 3 = 2.1666666666666665
  expect_equal(res$joint_stat, 13 / 6, tolerance = 1e-12)
})

# ============================================================================
# VM-7.5.2: F-test p-value verification
# Input: F = 13/6, df1 = 3, df2 = 48
# Expected: pf(13/6, 3, 48, lower.tail=FALSE) = 0.104158973181339
# R reference value computed independently
# ============================================================================
test_that("VM-7.5.2: F-test p-value verified against R pf() reference", {
  pe <- make_vm_effects(
    t_stats = c(2.0, 1.5, 0.5),
    n_treated = rep(25L, 3), n_control = rep(25L, 3)
  )
  res <- test_parallel_trends_joint(pe, test_type = "f")
  # df2 = max(as.integer(mean(25,25,25) + mean(25,25,25) - 2), 1) = 48
  expect_equal(res$joint_df2, 48L)
  # R reference: pf(13/6, 3, 48, lower.tail=FALSE) = 0.104158973181339
  expected_p <- 0.104158973181339
  expect_equal(res$joint_pvalue, expected_p, tolerance = 1e-12)
})

# ============================================================================
# VM-7.5.3: Wald statistic verification
# Input: t_stats = [2.0, 1.5, 0.5], K = 3
# Expected: W = 4 + 2.25 + 0.25 = 6.5 (exact)
# vibe-math result: 6.5
# ============================================================================
test_that("VM-7.5.3: Wald statistic = sum(t²) verified by vibe-math", {
  pe <- make_vm_effects(
    t_stats = c(2.0, 1.5, 0.5),
    n_treated = rep(25L, 3), n_control = rep(25L, 3)
  )
  res <- test_parallel_trends_joint(pe, test_type = "wald")
  # vibe-math: 2.0^2 + 1.5^2 + 0.5^2 = 6.5
  expect_equal(res$joint_stat, 6.5, tolerance = 1e-12)
})

# ============================================================================
# VM-7.5.4: Wald p-value verification
# Input: W = 6.5, df = 3
# Expected: pchisq(6.5, df=3, lower.tail=FALSE) = 0.089662503988168
# R reference value computed independently
# ============================================================================
test_that("VM-7.5.4: Wald p-value verified against R pchisq() reference", {
  pe <- make_vm_effects(
    t_stats = c(2.0, 1.5, 0.5),
    n_treated = rep(25L, 3), n_control = rep(25L, 3)
  )
  res <- test_parallel_trends_joint(pe, test_type = "wald")
  # R reference: pchisq(6.5, df=3, lower.tail=FALSE) = 0.089662503988168
  expected_p <- 0.089662503988168
  expect_equal(res$joint_pvalue, expected_p, tolerance = 1e-12)
})

# ============================================================================
# VM-7.5.5: df2 calculation verification
# Input: n_treated = [10, 12, 8], n_control = [20, 18, 22]
# Expected: avg_n = mean(10,12,8) + mean(20,18,22) = 10 + 20 = 30
#           df2 = max(as.integer(30 - 2), 1) = 28
# vibe-math result: 28.0
# ============================================================================
test_that("VM-7.5.5: df2 calculation verified by vibe-math", {
  pe <- make_vm_effects(
    t_stats = c(2.0, 1.5, 0.5),
    n_treated = c(10L, 12L, 8L),
    n_control = c(20L, 18L, 22L)
  )
  res <- test_parallel_trends_joint(pe, test_type = "f")
  # vibe-math: (10+12+8)/3 + (20+18+22)/3 - 2 = 28
  expect_equal(res$joint_df2, 28L)
})
