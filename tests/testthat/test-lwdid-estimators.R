# ============================================================================
# test-lwdid-estimators.R — E2E tests for lwdid() with all four estimators
# Story E6-05, Task E6-05.5: TC-6.5.30 to TC-6.5.36, TC-6.5.47
# ============================================================================

# --- Setup ---
if (!exists(".lwdid_env") ||
    is.null(.lwdid_env$warning_registry)) {
  .lwdid_env$warning_registry <- new_warning_registry()
}

# --- Helper: Common Timing data with controls ---
# Returns data with d (treatment indicator) and post (post-period indicator)
# for Common Timing mode, plus gvar for staggered-compatible calls.
make_ct_est_data <- function(n = 300, seed = 123) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  d <- rbinom(n, 1, plogis(0.3 * x1 - 0.2 * x2))
  ids <- seq_len(n)
  df <- data.frame(
    id = rep(ids, each = 2),
    year = rep(c(2000, 2001), n),
    d = rep(d, each = 2),
    post = rep(c(0L, 1L), n),
    gvar = ifelse(rep(d, each = 2) == 1, 2001, 0),
    x1 = rep(x1, each = 2),
    x2 = rep(x2, each = 2)
  )
  df$y <- 1 + 0.5 * df$x1 - 0.3 * df$x2 +
    0.5 * (df$year == 2001) +
    2.0 * (df$d == 1 & df$post == 1) +
    rnorm(nrow(df), sd = 0.5)
  df
}

# --- Helper: Staggered data with controls ---
make_stag_est_data <- function(n_per = 100, seed = 456) {
  set.seed(seed)
  n_total <- 3 * n_per
  ids <- seq_len(n_total)
  groups <- rep(c(0, 2002, 2003), each = n_per)
  x1 <- rnorm(n_total)
  x2 <- rnorm(n_total)
  df <- expand.grid(id = ids, year = 2000:2004)
  df <- df[order(df$id, df$year), ]
  df$gvar <- groups[df$id]
  df$x1 <- x1[df$id]
  df$x2 <- x2[df$id]
  df$y <- 1 + 0.5 * df$x1 - 0.3 * df$x2 +
    0.3 * (df$year - 2000) +
    2.0 * (df$gvar > 0 & df$year >= df$gvar) +
    rnorm(nrow(df), sd = 0.5)
  df
}


# ============================================================================
# TC-6.5.30: Common Timing + IPW
# ============================================================================
test_that("TC-6.5.30: Common Timing + IPW end-to-end", {
  df <- make_ct_est_data()
  res <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", controls = c("x1", "x2"),
          estimator = "ipw")
  ))
  expect_true(is.numeric(res$att))
  expect_true(is.finite(res$att))
  expect_true(res$se_att > 0)
  expect_true(res$ci_lower < res$ci_upper)
  expect_equal(res$estimator, "ipw")
})

# ============================================================================
# TC-6.5.31: Staggered + IPWRA
# ============================================================================
test_that("TC-6.5.31: Staggered + IPWRA end-to-end", {
  df <- make_stag_est_data()
  res <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          gvar = "gvar", controls = c("x1", "x2"),
          estimator = "ipwra")
  ))
  expect_true(res$is_staggered)
  expect_equal(res$estimator, "ipwra")
  # Should have group-time ATTs
  expect_true(is.data.frame(res$att_by_cohort_time))
  expect_true(nrow(res$att_by_cohort_time) > 0)
})

# ============================================================================
# TC-6.5.32: PSM + caliper
# ============================================================================
test_that("TC-6.5.32: Common Timing + PSM + caliper", {
  df <- make_ct_est_data()
  res <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", controls = c("x1", "x2"),
          estimator = "psm", caliper = 0.2)
  ))
  expect_true(is.numeric(res$att))
  expect_true(is.finite(res$att))
  expect_true(res$se_att > 0)
  expect_equal(res$estimator, "psm")
})


# ============================================================================
# TC-6.5.33: Staggered + IPW per (g,r)
# ============================================================================
test_that("TC-6.5.33: Staggered + IPW per (g,r)", {
  df <- make_stag_est_data()
  res <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          gvar = "gvar", controls = c("x1", "x2"),
          estimator = "ipw")
  ))
  expect_true(res$is_staggered)
  expect_equal(res$estimator, "ipw")
  expect_true(is.data.frame(res$att_by_cohort_time))
  expect_true(nrow(res$att_by_cohort_time) > 0)
})

# ============================================================================
# TC-6.5.34: Four estimators comparison on same data
# ============================================================================
test_that("TC-6.5.34: Four estimators ATT comparison", {
  df <- make_ct_est_data(n = 400, seed = 789)
  estimators <- c("ra", "ipw", "ipwra", "psm")
  atts <- ses <- numeric(4)
  for (i in seq_along(estimators)) {
    res <- suppressWarnings(suppressMessages(
      lwdid(data = df, y = "y", ivar = "id", tvar = "year",
            d = "d", post = "post", controls = c("x1", "x2"),
            estimator = estimators[i])
    ))
    atts[i] <- res$att
    ses[i] <- res$se_att
  }
  # All ATTs should be finite
  expect_true(all(is.finite(atts)))
  expect_true(all(ses > 0))
  # ATT range should be reasonable: range < 2 * max(SE)
  att_range <- max(atts) - min(atts)
  expect_true(att_range < 2 * max(ses),
    info = sprintf("ATT range=%.4f, 2*max(SE)=%.4f", att_range, 2 * max(ses)))
})

# ============================================================================
# TC-6.5.35: Staggered + IPW [R enhancement]
# ============================================================================
test_that("TC-6.5.35: Staggered + IPW works (R enhancement)", {
  df <- make_stag_est_data()
  # This should NOT error — R supports staggered+IPW unlike Python
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          gvar = "gvar", controls = c("x1", "x2"),
          estimator = "ipw")
  )))
})

# ============================================================================
# TC-6.5.36: boot_reps passthrough
# ============================================================================
test_that("TC-6.5.36: boot_reps passthrough with bootstrap SE", {
  df <- make_ct_est_data(n = 200, seed = 111)
  # Use small boot_reps for speed
  res <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", controls = c("x1", "x2"),
          estimator = "ipw", se_method = "bootstrap",
          boot_reps = 50, seed = 42)
  ))
  expect_true(is.numeric(res$att))
  expect_true(res$se_att > 0)
  expect_equal(res$estimator, "ipw")
})


# ============================================================================
# TC-6.5.47: RA regression compatibility
# estimator="ra" (explicit) should give identical ATT/SE/CI as default
# ============================================================================
test_that("TC-6.5.47: estimator='ra' regression compatibility", {
  df <- make_ct_est_data(n = 300, seed = 999)
  # Default call (estimator="ra" is the default)
  res_default <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", controls = c("x1", "x2"))
  ))
  # Explicit estimator="ra"
  res_ra <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", controls = c("x1", "x2"),
          estimator = "ra")
  ))
  # ATT, SE, CI must be exactly identical (difference = 0)
  expect_identical(res_default$att, res_ra$att)
  expect_identical(res_default$se_att, res_ra$se_att)
  expect_identical(res_default$ci_lower, res_ra$ci_lower)
  expect_identical(res_default$ci_upper, res_ra$ci_upper)
  expect_equal(res_ra$estimator, "ra")
})

# ============================================================================
# Numerical reasonableness: all four estimators ATT near true effect ~2.0
# ============================================================================
test_that("Numerical reasonableness: ATTs near true effect 2.0", {
  df <- make_ct_est_data(n = 500, seed = 2024)
  estimators <- c("ra", "ipw", "ipwra", "psm")
  atts <- numeric(4)
  for (i in seq_along(estimators)) {
    res <- suppressWarnings(suppressMessages(
      lwdid(data = df, y = "y", ivar = "id", tvar = "year",
            d = "d", post = "post", controls = c("x1", "x2"),
            estimator = estimators[i])
    ))
    atts[i] <- res$att
  }
  names(atts) <- estimators
  # True DGP effect is 2.0; each ATT should be within reasonable range
  # With n=500 and sd=0.5, ATTs should be within ~1.0 of true effect
  for (est in estimators) {
    expect_true(
      abs(atts[est] - 2.0) < 1.5,
      info = sprintf("estimator=%s: ATT=%.4f, expected near 2.0", est, atts[est])
    )
  }
  # Coefficient of variation of ATTs should be small
  att_range <- max(atts) - min(atts)
  att_mean <- mean(atts)
  expect_true(
    att_range / abs(att_mean) < 0.5,
    info = sprintf("ATT range/mean = %.4f/%.4f = %.4f",
                   att_range, att_mean, att_range / abs(att_mean))
  )
})
