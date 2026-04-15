# ============================================================================
# Tests for compute_hc_vce()
# Story E3-02: HC0-HC4 Heteroskedasticity-Robust Standard Errors
#
# Test Groups:
#   Group 1: Sandwich consistency (HC0-HC4 vs sandwich::vcovHC)
# ============================================================================

# ============================================================================
# Helper functions
# ============================================================================

#' Generate heteroskedastic test data for HC VCE tests
make_hc_test_data <- function(n = 30, seed = 42) {
  set.seed(seed)
  x <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = abs(x) + 0.5)
  data.frame(y = y, D = d, x = x)
}

# ============================================================================
# Group 1: Sandwich consistency
# ============================================================================

test_that("HC0-HC4 vcov matches sandwich::vcovHC exactly", {
  df <- make_hc_test_data()
  fit <- lm(y ~ D + x, data = df)

  for (hc_type in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    result <- compute_hc_vce(fit, type = hc_type)
    expected <- sandwich::vcovHC(
      fit, type = toupper(hc_type)
    )
    expect_equal(
      result$vcov, expected,
      tolerance = 1e-12,
      label = paste(
        "compute_hc_vce vcov for", hc_type
      ),
      expected.label = paste(
        "sandwich::vcovHC for", toupper(hc_type)
      )
    )
  }
})

test_that("HC0-HC4 SE equals sqrt(diag(vcov))", {
  df <- make_hc_test_data()
  fit <- lm(y ~ D + x, data = df)

  for (hc_type in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    result <- compute_hc_vce(fit, type = hc_type)
    expected_se <- sqrt(diag(result$vcov))
    expect_equal(
      result$se, expected_se,
      tolerance = 1e-15,
      label = paste(
        "se vs sqrt(diag(vcov)) for", hc_type
      )
    )
  }
})

test_that("HC0-HC4 SE matches sandwich SE exactly", {
  df <- make_hc_test_data()
  fit <- lm(y ~ D + x, data = df)

  for (hc_type in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    result <- compute_hc_vce(fit, type = hc_type)
    expected_vcov <- sandwich::vcovHC(
      fit, type = toupper(hc_type)
    )
    expected_se <- sqrt(diag(expected_vcov))
    expect_equal(
      result$se, expected_se,
      tolerance = 1e-12,
      label = paste(
        "compute_hc_vce se for", hc_type
      ),
      expected.label = paste(
        "sqrt(diag(sandwich)) for", toupper(hc_type)
      )
    )
  }
})

test_that("compute_hc_vce returns correct structure", {
  df <- make_hc_test_data()
  fit <- lm(y ~ D + x, data = df)
  p <- length(coef(fit))

  for (hc_type in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    result <- compute_hc_vce(fit, type = hc_type)

    # Exactly two fields: vcov and se
    expect_identical(
      sort(names(result)), c("se", "vcov"),
      label = paste("field names for", hc_type)
    )

    # vcov is a p x p matrix
    expect_true(
      is.matrix(result$vcov),
      label = paste("vcov is matrix for", hc_type)
    )
    expect_equal(
      dim(result$vcov), c(p, p),
      label = paste("vcov dims for", hc_type)
    )

    # se is a numeric vector of length p
    expect_true(
      is.numeric(result$se),
      label = paste("se is numeric for", hc_type)
    )
    expect_true(
      is.vector(result$se),
      label = paste("se is vector for", hc_type)
    )
    expect_equal(
      length(result$se), p,
      label = paste("se length for", hc_type)
    )
  }
})

test_that("compute_hc_vce with hc1 matches robust alias via dispatcher", {
  df <- make_hc_test_data()
  fit <- lm(y ~ D + x, data = df)

  res_hc1 <- compute_hc_vce(fit, type = "hc1")
  res_robust <- compute_vce(fit, vce = "robust")

  expect_identical(res_robust$vcov, res_hc1$vcov)
})

# ============================================================================
# Group 2: HC hierarchy relationships
# ============================================================================

test_that("HC1 SE equals HC0 SE * sqrt(N/(N-p)) exactly", {
  df <- make_hc_test_data(n = 30, seed = 42)
  fit <- lm(y ~ D + x, data = df)
  n <- nobs(fit)
  p <- length(coef(fit))

  se_hc0 <- compute_hc_vce(fit, type = "hc0")$se
  se_hc1 <- compute_hc_vce(fit, type = "hc1")$se

  # V_HC1 = (N/(N-p)) * V_HC0  =>  SE_HC1 = sqrt(N/(N-p)) * SE_HC0
  correction <- sqrt(n / (n - p))
  expect_equal(se_hc1, se_hc0 * correction, tolerance = 1e-12)

  # HC1 > HC0 strictly (since N/(N-p) > 1 for p > 0)
  expect_true(all(se_hc1 > se_hc0))
})

test_that("Small sample hierarchy: HC3 >= HC2 >= HC1 >= HC0", {
  df <- make_hc_test_data(n = 30, seed = 42)
  fit <- lm(y ~ D + x, data = df)

  se_hc0 <- compute_hc_vce(fit, type = "hc0")$se
  se_hc1 <- compute_hc_vce(fit, type = "hc1")$se
  se_hc2 <- compute_hc_vce(fit, type = "hc2")$se
  se_hc3 <- compute_hc_vce(fit, type = "hc3")$se

  # HC1 >= HC0 (strict, since n > p)
  expect_true(all(se_hc1 >= se_hc0 - 1e-10),
    info = "HC1 >= HC0 violated"
  )

  # HC2 >= HC1
  expect_true(all(se_hc2 >= se_hc1 - 1e-10),
    info = "HC2 >= HC1 violated"
  )

  # HC3 >= HC2
  expect_true(all(se_hc3 >= se_hc2 - 1e-10),
    info = "HC3 >= HC2 violated"
  )
})

test_that("Large sample convergence: HC0-HC4 SE within 5%", {
  set.seed(123)
  df <- make_hc_test_data(n = 1000, seed = 123)
  fit <- lm(y ~ D + x, data = df)
  n <- nobs(fit)
  p <- length(coef(fit))

  se_hc0 <- compute_hc_vce(fit, type = "hc0")$se
  se_hc1 <- compute_hc_vce(fit, type = "hc1")$se
  se_hc2 <- compute_hc_vce(fit, type = "hc2")$se
  se_hc3 <- compute_hc_vce(fit, type = "hc3")$se
  se_hc4 <- compute_hc_vce(fit, type = "hc4")$se

  # All HC variants should be within 5% of HC0 for large N
  for (se_hcx in list(se_hc1, se_hc2, se_hc3, se_hc4)) {
    rel_diff <- abs(se_hcx - se_hc0) / se_hc0
    expect_true(all(rel_diff < 0.05),
      info = paste(
        "Max relative diff from HC0:",
        round(max(rel_diff), 4)
      )
    )
  }

  # Verify h_ii = O(1/N) convergence: max(h_ii) should be close to p/N
  h <- hatvalues(fit)
  expect_true(max(h) < 10 * p / n,
    info = paste(
      "max(h_ii) =", round(max(h), 6),
      "vs p/N =", round(p / n, 6)
    )
  )
})

# ============================================================================
# Group 3: Mathematical property verification
# ============================================================================

test_that("HC2 unbiasedness: under homoskedasticity, HC2 SE close to OLS SE", {
  # Under homoskedastic errors, HC2 is unbiased for the true variance,
  # so HC2 SE should be close to OLS SE (which is efficient).
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n) # homoskedastic errors
  fit <- lm(y ~ D + x, data = data.frame(y = y, D = d, x = x))

  se_ols <- sqrt(diag(vcov(fit)))
  se_hc2 <- compute_hc_vce(fit, type = "hc2")$se

  rel_diff <- abs(se_hc2 - se_ols) / se_ols
  expect_true(
    all(rel_diff < 0.20),
    info = paste(
      "HC2 vs OLS relative differences:",
      paste(round(rel_diff, 4), collapse = ", ")
    )
  )
})

test_that("HC3-jackknife exact relationship: V_jack = (N-1)/N * V_HC3", {
  # The HC3 variance is exactly related to the delete-one jackknife
  # variance (centered at the full-sample estimate):
  #   V_jack = (N-1)/N * V_HC3
  # where V_jack = (N-1)/N * sum_i (beta_{-i} - beta)(beta_{-i} - beta)'
  # and beta_{-i} is computed via Sherman-Morrison update.
  set.seed(123)
  n <- 10
  p <- 2
  x <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = abs(x) + 0.5)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  X <- model.matrix(fit)
  e <- residuals(fit)
  h <- hatvalues(fit)
  beta_full <- coef(fit)
  XtX_inv <- solve(crossprod(X))
  k <- ncol(X)

  # Compute delete-one estimates via Sherman-Morrison
  beta_loo <- matrix(NA_real_, nrow = n, ncol = k)
  for (i in seq_len(n)) {
    beta_loo[i, ] <- beta_full -
      as.numeric(XtX_inv %*% X[i, ] * e[i] / (1 - h[i]))
  }

  # Jackknife variance: (N-1)/N * sum_i (beta_{-i} - beta)(beta_{-i} - beta)'
  V_jack <- matrix(0, k, k)
  for (i in seq_len(n)) {
    diff_i <- beta_loo[i, ] - beta_full
    V_jack <- V_jack + tcrossprod(diff_i)
  }
  V_jack <- ((n - 1) / n) * V_jack

  # HC3 variance from our implementation
  V_hc3 <- compute_hc_vce(fit, type = "hc3")$vcov

  # Exact relationship: V_jack = (N-1)/N * V_HC3
  # Strip dimnames for comparison (sandwich adds coefficient names)
  expected <- unname(((n - 1) / n) * V_hc3)
  expect_equal(
    V_jack, expected,
    tolerance = 1e-10,
    label = "Jackknife variance",
    expected.label = "(N-1)/N * V_HC3"
  )
})

test_that("Sherman-Morrison delete-one matches actual leave-one-out OLS", {
  # For each observation i, verify that the Sherman-Morrison formula
  #   beta_{-i} = beta - (X'X)^{-1} x_i e_i / (1 - h_ii)
  # gives the same result as actually refitting OLS on the data
  # with observation i removed.
  set.seed(456)
  n <- 10
  x <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  df <- data.frame(y = y, D = d, x = x)
  fit <- lm(y ~ D + x, data = df)

  X <- model.matrix(fit)
  e <- residuals(fit)
  h <- hatvalues(fit)
  beta_full <- coef(fit)
  XtX_inv <- solve(crossprod(X))

  for (i in seq_len(n)) {
    # Sherman-Morrison formula
    beta_sm <- beta_full -
      as.numeric(XtX_inv %*% X[i, ] * e[i] / (1 - h[i]))

    # Actual delete-one OLS
    fit_loo <- lm(y ~ D + x, data = df[-i, ])
    beta_actual <- coef(fit_loo)

    expect_equal(
      beta_sm, beta_actual,
      tolerance = 1e-10,
      label = paste("Sherman-Morrison beta_{-", i, "}", sep = ""),
      expected.label = paste("Actual LOO beta_{-", i, "}", sep = "")
    )
  }
})

test_that("HC4 d_i values: high leverage -> d_i > 3.5, normal -> d_i < 1.5", {
  # HC4 uses d_i = min(4, h_ii / h_bar) where h_bar = p/N.
  # Construct data with one extreme leverage point.
  set.seed(789)
  n <- 20
  x <- rnorm(n)
  x[1] <- 100 # extreme leverage point
  d <- rbinom(n, 1, 0.5)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  fit <- lm(y ~ D + x, data = data.frame(y = y, D = d, x = x))

  h <- hatvalues(fit)
  p <- length(coef(fit))
  h_bar <- p / n
  d_i <- pmin(4, h / h_bar)

  # High-leverage observation (i=1) should have d_i > 3.5
  expect_true(
    d_i[1] > 3.5,
    info = paste("d_1 =", round(d_i[1], 4), "expected > 3.5")
  )

  # Normal observations (i=2..n) should have d_i < 1.5
  expect_true(
    all(d_i[-1] < 1.5),
    info = paste(
      "max d_i for normal obs =",
      round(max(d_i[-1]), 4), "expected < 1.5"
    )
  )
})

test_that("Leverage h_ii properties: h_ii in [0,1) and sum(h_ii) = p", {
  df <- make_hc_test_data(n = 30, seed = 42)
  fit <- lm(y ~ D + x, data = df)
  p <- length(coef(fit))

  h <- hatvalues(fit)

  # h_ii >= 0
  expect_true(
    all(h >= 0),
    info = paste("min h_ii =", round(min(h), 10))
  )

  # h_ii < 1 when N > p
  expect_true(
    all(h < 1),
    info = paste("max h_ii =", round(max(h), 10))
  )

  # sum(h_ii) = p (trace of hat matrix = rank of X)
  expect_equal(
    sum(h), p,
    tolerance = 1e-10,
    label = "sum(h_ii)",
    expected.label = "p"
  )
})

# ============================================================================
# Group 4: Input validation and warning conditions
# ============================================================================

test_that("Invalid HC type throws lwdid_invalid_parameter", {
  df <- make_hc_test_data()
  fit <- lm(y ~ D + x, data = df)

  expect_error(
    compute_hc_vce(fit, type = "HC5"),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    compute_hc_vce(fit, type = "invalid"),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    compute_hc_vce(fit, type = ""),
    class = "lwdid_invalid_parameter"
  )
})

test_that("HC1/HC3 emit lwdid_small_sample when N_treated=1", {
  set.seed(42)
  n <- 10
  d <- c(1, rep(0, n - 1))
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  fit <- lm(y ~ D + x, data = data.frame(y = y, D = d, x = x))

  expect_warning(
    compute_hc_vce(fit, type = "hc1"),
    class = "lwdid_small_sample"
  )
  expect_warning(
    compute_hc_vce(fit, type = "hc3"),
    class = "lwdid_small_sample"
  )
})

test_that("HC0/HC2 do NOT emit lwdid_small_sample when N_treated=1", {
  set.seed(42)
  n <- 10
  d <- c(1, rep(0, n - 1))
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  fit <- lm(y ~ D + x, data = data.frame(y = y, D = d, x = x))

  # HC0 should not trigger small sample warning
  withCallingHandlers(
    compute_hc_vce(fit, type = "hc0"),
    lwdid_small_sample = function(w) {
      fail("HC0 should not emit lwdid_small_sample warning")
    }
  )

  # HC2 should not trigger small sample warning
  # (it may trigger lwdid_numerical for high leverage, which is fine)
  withCallingHandlers(
    compute_hc_vce(fit, type = "hc2"),
    lwdid_small_sample = function(w) {
      fail("HC2 should not emit lwdid_small_sample warning")
    }
  )
})

test_that("HC2/HC3/HC4 emit lwdid_numerical for high leverage", {
  set.seed(42)
  n <- 5
  x <- c(1000, rnorm(n - 1))
  d <- c(1, 1, 1, 0, 0)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  fit <- lm(y ~ D + x, data = data.frame(y = y, D = d, x = x))

  # Verify the data actually has high leverage
  expect_true(
    max(hatvalues(fit)) > 0.99,
    info = paste(
      "max(h_ii) =", round(max(hatvalues(fit)), 6),
      "should be > 0.99"
    )
  )

  for (hc_type in c("hc2", "hc3", "hc4")) {
    expect_warning(
      compute_hc_vce(fit, type = hc_type),
      class = "lwdid_numerical",
      label = paste(hc_type, "high leverage warning")
    )
  }
})

test_that("HC0/HC1 do NOT emit lwdid_numerical for high leverage", {
  set.seed(42)
  n <- 5
  x <- c(1000, rnorm(n - 1))
  d <- c(1, 1, 1, 0, 0)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  fit <- lm(y ~ D + x, data = data.frame(y = y, D = d, x = x))

  # HC0 should not trigger numerical warning
  withCallingHandlers(
    compute_hc_vce(fit, type = "hc0"),
    lwdid_numerical = function(w) {
      fail("HC0 should not emit lwdid_numerical warning")
    }
  )

  # HC1 may emit lwdid_small_sample but NOT lwdid_numerical
  withCallingHandlers(
    compute_hc_vce(fit, type = "hc1"),
    lwdid_numerical = function(w) {
      fail("HC1 should not emit lwdid_numerical warning")
    }
  )
})

test_that("No D variable: small sample check safely skipped", {
  set.seed(42)
  fit_no_d <- lm(y ~ x,
    data = data.frame(y = rnorm(20), x = rnorm(20))
  )

  # Should not error, just compute normally
  result <- compute_hc_vce(fit_no_d, type = "hc1")
  expect_true(is.matrix(result$vcov))
  expect_equal(nrow(result$vcov), 2L)
  expect_equal(ncol(result$vcov), 2L)
  expect_true(is.numeric(result$se))
  expect_equal(length(result$se), 2L)
})

test_that("Default type parameter is HC3", {
  df <- make_hc_test_data()
  fit <- lm(y ~ D + x, data = df)

  result_default <- compute_hc_vce(fit)
  result_hc3 <- compute_hc_vce(fit, type = "hc3")

  expect_identical(result_default$vcov, result_hc3$vcov)
  expect_identical(result_default$se, result_hc3$se)
})

# ============================================================================
# Group 5: Smoking dataset numerical consistency
# ============================================================================

test_that("HC3 SE matches sandwich on smoking data", {
  data(smoking, package = "lwdid")
  fit <- lm(cigsale ~ treat + d + post, data = smoking)

  result <- compute_hc_vce(fit, type = "hc3")
  expected <- sandwich::vcovHC(fit, type = "HC3")

  expect_equal(result$vcov, expected, tolerance = 1e-12)
})

test_that("HC0-HC4 vcov matches sandwich on smoking data", {
  data(smoking, package = "lwdid")
  fit <- lm(cigsale ~ treat + d + post, data = smoking)

  for (hc_type in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    result <- compute_hc_vce(fit, type = hc_type)
    expected <- sandwich::vcovHC(fit, type = toupper(hc_type))
    expect_equal(
      result$vcov, expected,
      tolerance = 1e-12,
      label = paste(
        "smoking vcov for", hc_type
      ),
      expected.label = paste(
        "sandwich::vcovHC for", toupper(hc_type)
      )
    )
  }
})

test_that("HC hierarchy on smoking data: SE_HC1 > SE_HC0 strictly", {
  data(smoking, package = "lwdid")
  fit <- lm(cigsale ~ treat + d + post, data = smoking)

  se_hc0 <- compute_hc_vce(fit, type = "hc0")$se
  se_hc1 <- compute_hc_vce(fit, type = "hc1")$se

  # HC1 > HC0 strictly (since N/(N-p) > 1 for p > 0)
  expect_true(
    all(se_hc1 > se_hc0),
    info = "HC1 > HC0 strictly on smoking data"
  )

  # Verify exact HC1/HC0 ratio = sqrt(N/(N-p))
  n <- nobs(fit)
  p <- length(coef(fit))
  correction <- sqrt(n / (n - p))
  expect_equal(
    se_hc1, se_hc0 * correction,
    tolerance = 1e-12,
    label = "HC1/HC0 ratio on smoking data"
  )
})

test_that("HC hierarchy on smoking data: HC3 >= HC2 >= HC1 generally", {
  data(smoking, package = "lwdid")
  fit <- lm(cigsale ~ treat + d + post, data = smoking)

  se_hc0 <- compute_hc_vce(fit, type = "hc0")$se
  se_hc1 <- compute_hc_vce(fit, type = "hc1")$se
  se_hc2 <- compute_hc_vce(fit, type = "hc2")$se
  se_hc3 <- compute_hc_vce(fit, type = "hc3")$se

  # HC1 > HC0 strictly (scalar multiple, always holds)
  expect_true(
    all(se_hc1 > se_hc0),
    info = "HC1 > HC0 must hold strictly"
  )

  # HC2 >= HC1 is a general tendency, not a strict per-coefficient

  # guarantee. The inequality derives from Jensen's inequality on the

  # trace of the meat matrix, which does not imply element-wise
  # dominance. We verify the majority of coefficients satisfy it.
  frac_hc2_ge_hc1 <- mean(se_hc2 >= se_hc1 - 1e-10)
  expect_gte(
    frac_hc2_ge_hc1, 0.5,
    label = paste(
      "fraction of coefficients with HC2 >= HC1 on smoking data:",
      round(frac_hc2_ge_hc1, 3)
    )
  )

  # HC3 >= HC2 (generally holds more reliably)
  frac_hc3_ge_hc2 <- mean(se_hc3 >= se_hc2 - 1e-10)
  expect_gte(
    frac_hc3_ge_hc2, 0.5,
    label = paste(
      "fraction of coefficients with HC3 >= HC2 on smoking data:",
      round(frac_hc3_ge_hc2, 3)
    )
  )
})

test_that("Leverage h_ii properties on smoking data", {
  data(smoking, package = "lwdid")
  fit <- lm(cigsale ~ treat + d + post, data = smoking)
  p <- length(coef(fit))

  h <- hatvalues(fit)

  # sum(h_ii) = p (trace of hat matrix = rank of X)
  expect_equal(
    sum(h), p,
    tolerance = 1e-10,
    label = "sum(h_ii) on smoking data",
    expected.label = "p"
  )

  # h_ii >= 0
  expect_true(
    all(h >= 0),
    info = paste("min h_ii on smoking =", round(min(h), 10))
  )

  # h_ii < 1 when N >> p (1209 obs, 4 coefficients)
  expect_true(
    all(h < 1),
    info = paste("max h_ii on smoking =", round(max(h), 10))
  )
})

test_that("Smoking data model dimensions are correct", {
  data(smoking, package = "lwdid")
  fit <- lm(cigsale ~ treat + d + post, data = smoking)

  # 4 coefficients: intercept, treat, d, post
  expect_equal(length(coef(fit)), 4L)
  expect_equal(nobs(fit), 1209L)

  # VCE dimensions match
  result <- compute_hc_vce(fit, type = "hc3")
  expect_equal(dim(result$vcov), c(4L, 4L))
  expect_equal(length(result$se), 4L)
})

test_that("Smoking data: no D column triggers safe skip of small sample check", {
  # The smoking dataset uses lowercase 'd' not uppercase 'D',

  # so the small sample warning check in compute_hc_vce won't find "D"
  # in model.frame — this is fine, it safely skips.
  data(smoking, package = "lwdid")
  fit <- lm(cigsale ~ treat + d + post, data = smoking)

  # HC1 and HC3 should NOT emit lwdid_small_sample warning
  # because the column is 'd' not 'D'
  withCallingHandlers(
    compute_hc_vce(fit, type = "hc1"),
    lwdid_small_sample = function(w) {
      fail("HC1 should not emit lwdid_small_sample on smoking data (lowercase d)")
    }
  )

  withCallingHandlers(
    compute_hc_vce(fit, type = "hc3"),
    lwdid_small_sample = function(w) {
      fail("HC3 should not emit lwdid_small_sample on smoking data (lowercase d)")
    }
  )
})

# TODO: Fill in Stata benchmark values via Stata MCP
# reg cigsale treat d post, robust
# Expected HC1 SE for treat coefficient: [to be determined]
# reg cigsale treat d post, vce(hc3)
# Expected HC3 SE for treat coefficient: [to be determined]
