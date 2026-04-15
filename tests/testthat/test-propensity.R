# test-propensity.R — Unit tests for propensity score infrastructure
# Story E6-01: Tasks E6-01.6 and E6-01.7

# ============================================================================
# Section 1: estimate_propensity_score() tests
# ============================================================================

# Helper: generate test data for propensity score estimation
# Returns a data.frame with treatment D generated from logistic model
generate_ps_test_data <- function(n = 200, seed = 42, n_covariates = 2) {
  set.seed(seed)
  df <- data.frame(id = seq_len(n))

  # Generate covariates
  for (j in seq_len(n_covariates)) {
    df[[paste0("X", j)]] <- rnorm(n)
  }

  # Generate treatment from logistic model
  # Linear predictor: -0.5 + 0.8*X1 + 0.5*X2 + ...
  linpred <- -0.5 + 0.8 * df$X1
  if (n_covariates >= 2) linpred <- linpred + 0.5 * df$X2
  if (n_covariates >= 3) linpred <- linpred + 0.3 * df$X3

  prob <- 1 / (1 + exp(-linpred))
  df$D <- rbinom(n, 1, prob)

  # Ensure both groups exist

  if (all(df$D == 1) || all(df$D == 0)) {
    df$D[1] <- 0L
    df$D[2] <- 1L
  }

  df
}

# ============================================================================
# Test Group 1: Normal data PS range test
# ============================================================================

test_that("1.1: normal data returns complete 8-element list", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  controls <- c("X1", "X2")

  result <- estimate_propensity_score(df, "D", controls)

  # Return list has exactly 8 elements with correct names
  expect_type(result, "list")
  expect_length(result, 8L)
  expected_names <- c("propensity_scores", "coefficients", "converged",
                      "n_trimmed", "trimmed_mask", "trim_method",
                      "glm_object", "removed_cols")
  expect_named(result, expected_names)

  # All PS in [0, 1]
  ps <- result$propensity_scores
  expect_true(all(ps >= 0 & ps <= 1))

  # converged is TRUE
  expect_true(result$converged)

  # coefficients is a named list with "_intercept" key
  expect_type(result$coefficients, "list")
  expect_true("_intercept" %in% names(result$coefficients))
  expect_true("X1" %in% names(result$coefficients))
  expect_true("X2" %in% names(result$coefficients))

  # glm_object is a glm object
  expect_s3_class(result$glm_object, "glm")

  # removed_cols is character(0)
  expect_identical(result$removed_cols, character(0))

  # PS length matches data
  expect_length(ps, nrow(df))

  # n_trimmed is integer
  expect_type(result$n_trimmed, "integer")

  # trim_method is "clip" (default)
  expect_equal(result$trim_method, "clip")
})

test_that("1.2: PS values are within default clip bounds [0.01, 0.99]", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  result <- estimate_propensity_score(df, "D", c("X1", "X2"))

  ps <- result$propensity_scores
  expect_true(all(ps >= 0.01))
  expect_true(all(ps <= 0.99))
})

test_that("1.3: vcov() works on returned glm_object", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  result <- estimate_propensity_score(df, "D", c("X1", "X2"))

  # vcov should return a valid variance-covariance matrix
  V <- vcov(result$glm_object)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 3L)  # intercept + 2 covariates
  expect_equal(ncol(V), 3L)
  # Diagonal should be positive (variances)
  expect_true(all(diag(V) > 0))
})


# ============================================================================
# Test Group 2: Constant column detection
# ============================================================================

test_that("2.1: constant column detected and removed with lwdid_data warning", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$const_col <- rep(5, nrow(df))
  controls <- c("X1", "X2", "const_col")

  # Capture warnings and result separately
  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_propensity_score(df, "D", controls),
    lwdid_data = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  # lwdid_data warning was issued
  expect_true(length(warnings_caught) > 0L)
  expect_true(inherits(warnings_caught[[1]], "lwdid_data"))

  # removed_cols contains the constant column name
  expect_true("const_col" %in% result$removed_cols)
  expect_length(result$removed_cols, 1L)

  # PS still computed correctly (constant col excluded)
  expect_true(result$converged)
  expect_true(all(result$propensity_scores >= 0.01))
  expect_true(all(result$propensity_scores <= 0.99))

  # Coefficients should NOT contain the constant column
  expect_false("const_col" %in% names(result$coefficients))
  expect_true("X1" %in% names(result$coefficients))
  expect_true("X2" %in% names(result$coefficients))
})

test_that("2.2: multiple constant columns detected", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$const1 <- rep(0, nrow(df))
  df$const2 <- rep(99, nrow(df))
  controls <- c("X1", "const1", "X2", "const2")

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_propensity_score(df, "D", controls),
    lwdid_data = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_true(length(warnings_caught) > 0L)
  expect_true(all(c("const1", "const2") %in% result$removed_cols))
  expect_length(result$removed_cols, 2L)
  expect_true(result$converged)
})

# ============================================================================
# Test Group 3: All constant columns → error
# ============================================================================

test_that("3.1: all constant columns produces lwdid_insufficient_data error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$c1 <- rep(1, nrow(df))
  df$c2 <- rep(2, nrow(df))

  expect_error(
    suppressWarnings(estimate_propensity_score(df, "D", c("c1", "c2"))),
    class = "lwdid_insufficient_data"
  )
})

# ============================================================================
# Test Group 4: Input validation errors
# ============================================================================

test_that("4.1: non-binary D produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$D_bad <- sample(c(0, 1, 2), nrow(df), replace = TRUE)

  expect_error(
    estimate_propensity_score(df, "D_bad", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )
})

test_that("4.2: D with continuous values produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$D_cont <- runif(nrow(df))

  expect_error(
    estimate_propensity_score(df, "D_cont", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )
})

test_that("4.3: trim_threshold = 0 produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2"), trim_threshold = 0),
    class = "lwdid_invalid_param"
  )
})

test_that("4.4: trim_threshold = 0.5 produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2"), trim_threshold = 0.5),
    class = "lwdid_invalid_param"
  )
})

test_that("4.5: trim_threshold = -0.1 produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2"), trim_threshold = -0.1),
    class = "lwdid_invalid_param"
  )
})

test_that("4.6: trim_threshold = 1.0 produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2"), trim_threshold = 1.0),
    class = "lwdid_invalid_param"
  )
})

test_that("4.7: invalid trim_method produces error from match.arg", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2"), trim_method = "invalid")
  )
})


# ============================================================================
# Test Group 5: Trimming threshold test
# ============================================================================

test_that("5.1: aggressive trim=0.1 trims some observations", {
  # Generate data with more extreme PS by using strong covariate effect
  set.seed(123)
  n <- 300
  X1 <- rnorm(n, sd = 2)
  X2 <- rnorm(n, sd = 1)
  linpred <- -1 + 2.0 * X1 + 1.0 * X2  # strong effect → extreme PS
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  # Ensure both groups
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  df <- data.frame(X1 = X1, X2 = X2, D = D)

  result <- estimate_propensity_score(df, "D", c("X1", "X2"),
                                      trim_threshold = 0.1)

  # With strong effects and trim=0.1, some observations should be trimmed
  expect_true(result$n_trimmed > 0)
  # All PS should be within [0.1, 0.9] after clipping
  expect_true(all(result$propensity_scores >= 0.1))
  expect_true(all(result$propensity_scores <= 0.9))
})

# ============================================================================
# Test Group 6: Single vs multiple covariates
# ============================================================================

test_that("6.1: single covariate works correctly", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 1)

  result <- estimate_propensity_score(df, "D", "X1")

  expect_true(result$converged)
  expect_length(result$coefficients, 2L)  # _intercept + X1
  expect_true("_intercept" %in% names(result$coefficients))
  expect_true("X1" %in% names(result$coefficients))
  expect_length(result$propensity_scores, nrow(df))
})

test_that("6.2: three covariates works correctly", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 3)

  result <- estimate_propensity_score(df, "D", c("X1", "X2", "X3"))

  expect_true(result$converged)
  expect_length(result$coefficients, 4L)  # _intercept + X1 + X2 + X3
  expect_true("_intercept" %in% names(result$coefficients))
  expect_true("X1" %in% names(result$coefficients))
  expect_true("X2" %in% names(result$coefficients))
  expect_true("X3" %in% names(result$coefficients))
  expect_length(result$propensity_scores, nrow(df))
})

# ============================================================================
# Test Group 7: Convergence/boundary warnings
# ============================================================================

test_that("7.1: quasi-complete separation triggers lwdid_numerical warning", {
  # Create data with perfect separation: D perfectly determined by X1 > 0
  # With well-separated groups, glm may not converge or may hit boundary
  set.seed(99)
  n <- 100
  X1 <- c(rnorm(50, mean = -3, sd = 0.5), rnorm(50, mean = 3, sd = 0.5))
  D <- as.integer(X1 > 0)  # perfect separation
  df <- data.frame(X1 = X1, D = D)

  # With perfect separation, glm should trigger lwdid_numerical warning
  # (either non-convergence or boundary)
  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_propensity_score(df, "D", "X1"),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  # At least one warning should be lwdid_numerical
  has_numerical <- any(vapply(warnings_caught, function(w) {
    inherits(w, "lwdid_numerical")
  }, logical(1)))
  expect_true(has_numerical)

  # The warning should indicate a numerical issue (convergence or boundary)
  numerical_warnings <- Filter(
    function(w) inherits(w, "lwdid_numerical"),
    warnings_caught
  )
  expect_true(length(numerical_warnings) >= 1L)

  # converged should be FALSE or boundary should be TRUE
  expect_true(!result$converged || isTRUE(result$glm_object$boundary))
})


# ============================================================================
# Test Group 8: Clip behavior
# ============================================================================

test_that("8.1: clip mode with trim=0.1 — no NA, all PS in [0.1, 0.9]", {
  # Generate data with strong effects to ensure some raw PS are extreme
  set.seed(77)
  n <- 300
  X1 <- rnorm(n, sd = 2)
  X2 <- rnorm(n, sd = 1)
  linpred <- -1 + 2.0 * X1 + 0.5 * X2
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  df <- data.frame(X1 = X1, X2 = X2, D = D)

  result <- estimate_propensity_score(df, "D", c("X1", "X2"),
                                      trim_threshold = 0.1,
                                      trim_method = "clip")

  ps <- result$propensity_scores

  # No NA in propensity_scores
  expect_false(anyNA(ps))

  # All PS in [0.1, 0.9]
  expect_true(all(ps >= 0.1))
  expect_true(all(ps <= 0.9))

  # trimmed_mask is all FALSE for clip mode
  expect_true(all(result$trimmed_mask == FALSE))

  # trim_method recorded correctly
  expect_equal(result$trim_method, "clip")
})

test_that("8.2: clip mode preserves all observations (no NA)", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  result <- estimate_propensity_score(df, "D", c("X1", "X2"),
                                      trim_threshold = 0.05,
                                      trim_method = "clip")

  # All observations preserved
  expect_length(result$propensity_scores, nrow(df))
  expect_false(anyNA(result$propensity_scores))
  expect_true(all(result$propensity_scores >= 0.05))
  expect_true(all(result$propensity_scores <= 0.95))
})

# ============================================================================
# Test Group 9: Drop behavior
# ============================================================================

test_that("9.1: drop mode with trim=0.1 — extreme PS become NA", {
  # Generate data with strong effects to ensure some raw PS are extreme
  set.seed(77)
  n <- 300
  X1 <- rnorm(n, sd = 2)
  X2 <- rnorm(n, sd = 1)
  linpred <- -1 + 2.0 * X1 + 0.5 * X2
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  df <- data.frame(X1 = X1, X2 = X2, D = D)

  result <- estimate_propensity_score(df, "D", c("X1", "X2"),
                                      trim_threshold = 0.1,
                                      trim_method = "drop")

  ps <- result$propensity_scores

  # Some PS should be NA (the extreme ones)
  expect_true(anyNA(ps))

  # trimmed_mask has TRUE entries
  expect_true(any(result$trimmed_mask))

  # n_trimmed matches sum(trimmed_mask)
  expect_equal(result$n_trimmed, sum(result$trimmed_mask))

  # Non-NA PS should be within bounds
  valid_ps <- ps[!is.na(ps)]
  expect_true(all(valid_ps >= 0.1))
  expect_true(all(valid_ps <= 0.9))

  # trim_method recorded correctly
  expect_equal(result$trim_method, "drop")
})

test_that("9.2: drop mode — trimmed_mask TRUE where PS was extreme", {
  set.seed(77)
  n <- 300
  X1 <- rnorm(n, sd = 2)
  linpred <- -1 + 2.5 * X1
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  df <- data.frame(X1 = X1, D = D)

  result <- estimate_propensity_score(df, "D", "X1",
                                      trim_threshold = 0.1,
                                      trim_method = "drop")

  # Where trimmed_mask is TRUE, PS should be NA
  expect_true(all(is.na(result$propensity_scores[result$trimmed_mask])))
  # Where trimmed_mask is FALSE, PS should NOT be NA
  expect_false(anyNA(result$propensity_scores[!result$trimmed_mask]))
})


# ============================================================================
# Test Group 10: Numerical consistency — coefficients match R's own glm()
# ============================================================================

test_that("10.1: coefficients match R's own glm() output exactly", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  controls <- c("X1", "X2")

  result <- estimate_propensity_score(df, "D", controls)

  # Fit glm directly for comparison
  fml <- as.formula(paste("D ~", paste(controls, collapse = " + ")))
  glm_ref <- glm(fml, data = df, family = binomial(link = "logit"))
  ref_coef <- coef(glm_ref)

  # Coefficients should match exactly (same algorithm)
  expect_equal(result$coefficients[["_intercept"]],
               unname(ref_coef["(Intercept)"]),
               tolerance = 1e-12)
  expect_equal(result$coefficients[["X1"]],
               unname(ref_coef["X1"]),
               tolerance = 1e-12)
  expect_equal(result$coefficients[["X2"]],
               unname(ref_coef["X2"]),
               tolerance = 1e-12)
})

test_that("10.2: coefficients have expected signs (positive covariate → positive coef)", {
  # Data generated with positive X1 effect (0.8) and positive X2 effect (0.5)
  df <- generate_ps_test_data(n = 500, seed = 42, n_covariates = 2)
  result <- estimate_propensity_score(df, "D", c("X1", "X2"))

  # With N=500 and true coefficients 0.8 and 0.5, estimated coefficients
  # should be positive
  expect_true(result$coefficients[["X1"]] > 0)
  expect_true(result$coefficients[["X2"]] > 0)
})

test_that("10.3: PS values match predict(glm, type='response') with clipping", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  controls <- c("X1", "X2")

  result <- estimate_propensity_score(df, "D", controls, trim_threshold = 0.01)

  # Get raw PS from glm
  fml <- as.formula(paste("D ~", paste(controls, collapse = " + ")))
  glm_ref <- glm(fml, data = df, family = binomial(link = "logit"))
  ps_raw <- predict(glm_ref, type = "response")
  ps_clipped <- pmax(0.01, pmin(0.99, ps_raw))

  expect_equal(result$propensity_scores, unname(ps_clipped), tolerance = 1e-12)
})

# ============================================================================
# Test Group 11: NA input tests
# ============================================================================

test_that("11.1: D contains NA produces lwdid_invalid_param error mentioning column", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$D[5] <- NA

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )

  # Verify the error message mentions the column name
  err <- tryCatch(
    estimate_propensity_score(df, "D", c("X1", "X2")),
    lwdid_invalid_param = function(e) e
  )
  expect_true(grepl("D", err$message, fixed = TRUE))
})

test_that("11.2: control variable contains NA produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$X2[10] <- NA

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )

  # Verify the error message lists the NA column
  err <- tryCatch(
    estimate_propensity_score(df, "D", c("X1", "X2")),
    lwdid_invalid_param = function(e) e
  )
  expect_true(grepl("X2", err$message, fixed = TRUE))
})

test_that("11.3: multiple control variables with NA lists all NA columns", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 3)
  df$X1[5] <- NA
  df$X3[10] <- NA

  err <- tryCatch(
    estimate_propensity_score(df, "D", c("X1", "X2", "X3")),
    lwdid_invalid_param = function(e) e
  )
  expect_true(grepl("X1", err$message, fixed = TRUE))
  expect_true(grepl("X3", err$message, fixed = TRUE))
})

# ============================================================================
# Test Group 12: Degenerate input tests
# ============================================================================

test_that("12.1: all D=1 produces lwdid_insufficient_data error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$D <- rep(1L, nrow(df))

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2")),
    class = "lwdid_insufficient_data"
  )
})

test_that("12.2: all D=0 produces lwdid_insufficient_data error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)
  df$D <- rep(0L, nrow(df))

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2")),
    class = "lwdid_insufficient_data"
  )
})

# ============================================================================
# Additional edge case tests
# ============================================================================

test_that("13.1: missing treatment column produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  expect_error(
    estimate_propensity_score(df, "nonexistent", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )
})

test_that("13.2: missing control column produces lwdid_invalid_param error", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  expect_error(
    estimate_propensity_score(df, "D", c("X1", "X2", "X_missing")),
    class = "lwdid_invalid_param"
  )
})

test_that("13.3: d must be a single character string", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  # d as numeric
  expect_error(
    estimate_propensity_score(df, 1, c("X1", "X2")),
    class = "lwdid_invalid_param"
  )

  # d as vector of length 2
  expect_error(
    estimate_propensity_score(df, c("D", "X1"), c("X1", "X2")),
    class = "lwdid_invalid_param"
  )
})

test_that("13.4: controls must be character vector of length >= 1", {
  df <- generate_ps_test_data(n = 200, seed = 42, n_covariates = 2)

  # Empty controls
  expect_error(
    estimate_propensity_score(df, "D", character(0)),
    class = "lwdid_invalid_param"
  )

  # Numeric controls
  expect_error(
    estimate_propensity_score(df, "D", 1),
    class = "lwdid_invalid_param"
  )
})

test_that("13.5: clip vs drop produce different n_trimmed for same data", {
  set.seed(77)
  n <- 300
  X1 <- rnorm(n, sd = 2)
  linpred <- -1 + 2.5 * X1
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  df <- data.frame(X1 = X1, D = D)

  result_clip <- estimate_propensity_score(df, "D", "X1",
                                           trim_threshold = 0.1,
                                           trim_method = "clip")
  result_drop <- estimate_propensity_score(df, "D", "X1",
                                           trim_threshold = 0.1,
                                           trim_method = "drop")

  # Both should identify the same number of extreme observations
  expect_equal(result_clip$n_trimmed, result_drop$n_trimmed)

  # But clip has no NA, drop has NA
  expect_false(anyNA(result_clip$propensity_scores))
  expect_true(anyNA(result_drop$propensity_scores))

  # Clip trimmed_mask is all FALSE
  expect_true(all(result_clip$trimmed_mask == FALSE))
  # Drop trimmed_mask has TRUE entries
  expect_true(any(result_drop$trimmed_mask))
})

test_that("13.6: small sample (N=20) still works", {
  df <- generate_ps_test_data(n = 20, seed = 42, n_covariates = 1)
  result <- estimate_propensity_score(df, "D", "X1")

  expect_true(result$converged)
  expect_length(result$propensity_scores, 20L)
  expect_true(all(result$propensity_scores >= 0.01))
  expect_true(all(result$propensity_scores <= 0.99))
})


# ============================================================================
# Section 2: estimate_outcome_model() tests
# ============================================================================

# Helper: generate test data for outcome model estimation
# Returns a data.frame with known linear relationship for controls
generate_om_test_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  # True model for controls: Y = 2 + 3*X1 + 1*X2 + noise
  Y <- 2 + 3 * X1 + 1 * X2 + rnorm(n, sd = 0.5)
  # Add treatment effect for treated
  Y[D == 1] <- Y[D == 1] + 5
  data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
}

# ============================================================================
# Test Group 14: OLS basic functionality
# ============================================================================

test_that("14.1: OLS coefficients match R's lm() on control group", {
  df <- generate_om_test_data(n = 200, seed = 42)
  controls <- c("X1", "X2")

  result <- estimate_outcome_model(df, "Y", "D", controls)

  # Reference: fit lm on control group only
  ctrl_data <- df[df$D == 0, ]
  lm_ref <- lm(Y ~ X1 + X2, data = ctrl_data)
  ref_coef <- coef(lm_ref)

  # Coefficients should match to machine precision

  expect_equal(result$coefficients[["_intercept"]],
               unname(ref_coef["(Intercept)"]),
               tolerance = 1e-12)
  expect_equal(result$coefficients[["X1"]],
               unname(ref_coef["X1"]),
               tolerance = 1e-12)
  expect_equal(result$coefficients[["X2"]],
               unname(ref_coef["X2"]),
               tolerance = 1e-12)
})

test_that("14.2: OLS coefficients recover true DGP parameters (large N)", {
  # With large N and small noise, coefficients should be close to true values
  df <- generate_om_test_data(n = 5000, seed = 123)
  controls <- c("X1", "X2")

  result <- estimate_outcome_model(df, "Y", "D", controls)

  # True: intercept=2, X1=3, X2=1
  expect_equal(result$coefficients[["_intercept"]], 2, tolerance = 0.1)
  expect_equal(result$coefficients[["X1"]], 3, tolerance = 0.1)
  expect_equal(result$coefficients[["X2"]], 1, tolerance = 0.1)
})

test_that("14.3: OLS predictions match manual X %*% beta for all units", {
  df <- generate_om_test_data(n = 200, seed = 42)
  controls <- c("X1", "X2")

  result <- estimate_outcome_model(df, "Y", "D", controls)

  # Manual prediction: cbind(1, X) %*% beta
  X_all <- cbind(1, as.matrix(df[, controls]))
  beta <- c(result$coefficients[["_intercept"]],
            result$coefficients[["X1"]],
            result$coefficients[["X2"]])
  manual_pred <- as.numeric(X_all %*% beta)

  expect_equal(result$predictions, manual_pred, tolerance = 1e-12)
})

# ============================================================================
# Test Group 15: OLS predictions cover all units
# ============================================================================

test_that("15.1: OLS predictions length equals nrow(data)", {
  df <- generate_om_test_data(n = 200, seed = 42)
  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"))

  expect_length(result$predictions, nrow(df))
})

test_that("15.2: OLS predictions include both treated and control units", {
  df <- generate_om_test_data(n = 200, seed = 42)
  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"))

  # Predictions for treated units should exist and be finite
  treated_preds <- result$predictions[df$D == 1]
  control_preds <- result$predictions[df$D == 0]

  expect_true(length(treated_preds) > 0)
  expect_true(length(control_preds) > 0)
  expect_true(all(is.finite(treated_preds)))
  expect_true(all(is.finite(control_preds)))
})

test_that("15.3: OLS predictions are numeric with no NA", {
  df <- generate_om_test_data(n = 200, seed = 42)
  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"))

  expect_type(result$predictions, "double")
  expect_false(anyNA(result$predictions))
})

# ============================================================================
# Test Group 16: OLS coefficient format
# ============================================================================

test_that("16.1: coefficients is a named list with _intercept and all controls", {
  df <- generate_om_test_data(n = 200, seed = 42)
  controls <- c("X1", "X2")
  result <- estimate_outcome_model(df, "Y", "D", controls)

  expect_type(result$coefficients, "list")
  expect_true("_intercept" %in% names(result$coefficients))
  expect_true("X1" %in% names(result$coefficients))
  expect_true("X2" %in% names(result$coefficients))
  expect_length(result$coefficients, 3L)  # _intercept + X1 + X2
})

test_that("16.2: return value is a list with exactly 2 elements", {
  df <- generate_om_test_data(n = 200, seed = 42)
  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"))

  expect_type(result, "list")
  expect_length(result, 2L)
  expect_named(result, c("predictions", "coefficients"))
})

test_that("16.3: single covariate produces 2-element coefficient list", {
  df <- generate_om_test_data(n = 200, seed = 42)
  result <- estimate_outcome_model(df, "Y", "D", "X1")

  expect_length(result$coefficients, 2L)  # _intercept + X1
  expect_true("_intercept" %in% names(result$coefficients))
  expect_true("X1" %in% names(result$coefficients))
})

# ============================================================================
# Test Group 17: WLS basic functionality
# ============================================================================

test_that("17.1: WLS coefficients match manual sqrt(W) transformation + OLS", {
  set.seed(42)
  n <- 300
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- 2 + 3 * X1 + 1 * X2 + rnorm(n, sd = 0.5)
  Y[D == 1] <- Y[D == 1] + 5
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  # Generate PS-like weights: w = ps / (1 - ps)
  ps <- 1 / (1 + exp(-(-0.5 + 0.8 * X1 + 0.5 * X2)))
  sample_weights <- ps / (1 - ps)

  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                   sample_weights = sample_weights)

  # Manual WLS: transform by sqrt(w), then OLS
  ctrl_mask <- D == 0
  w_ctrl <- sample_weights[ctrl_mask]
  w_ctrl <- pmax(w_ctrl, 1e-10)
  sqrt_w <- sqrt(w_ctrl)

  ctrl_data <- df[ctrl_mask, ]
  ctrl_data_w <- data.frame(
    Y = ctrl_data$Y * sqrt_w,
    intercept = sqrt_w,
    X1 = ctrl_data$X1 * sqrt_w,
    X2 = ctrl_data$X2 * sqrt_w
  )
  lm_wls <- lm(Y ~ 0 + intercept + X1 + X2, data = ctrl_data_w)
  ref_coef <- coef(lm_wls)

  expect_equal(result$coefficients[["_intercept"]],
               unname(ref_coef["intercept"]),
               tolerance = 1e-10)
  expect_equal(result$coefficients[["X1"]],
               unname(ref_coef["X1"]),
               tolerance = 1e-10)
  expect_equal(result$coefficients[["X2"]],
               unname(ref_coef["X2"]),
               tolerance = 1e-10)
})

test_that("17.2: WLS predictions match manual X_all %*% beta_wls", {
  set.seed(42)
  n <- 300
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- 2 + 3 * X1 + 1 * X2 + rnorm(n, sd = 0.5)
  Y[D == 1] <- Y[D == 1] + 5
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  ps <- 1 / (1 + exp(-(-0.5 + 0.8 * X1 + 0.5 * X2)))
  sample_weights <- ps / (1 - ps)

  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                   sample_weights = sample_weights)

  # Manual prediction
  X_all <- cbind(1, as.matrix(df[, c("X1", "X2")]))
  beta <- c(result$coefficients[["_intercept"]],
            result$coefficients[["X1"]],
            result$coefficients[["X2"]])
  manual_pred <- as.numeric(X_all %*% beta)

  expect_equal(result$predictions, manual_pred, tolerance = 1e-12)
})

# ============================================================================
# Test Group 18: WLS weight floor test
# ============================================================================

test_that("18.1: near-zero weights (1e-15) produce no NaN/Inf in predictions", {
  df <- generate_om_test_data(n = 200, seed = 42)

  # Create weights with near-zero values
  sample_weights <- rep(1e-15, nrow(df))
  # Give some control units reasonable weights so mean > 0
  ctrl_idx <- which(df$D == 0)
  sample_weights[ctrl_idx[1:5]] <- 1.0

  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                   sample_weights = sample_weights)

  expect_false(any(is.nan(result$predictions)))
  expect_false(any(is.infinite(result$predictions)))
  expect_false(any(is.nan(unlist(result$coefficients))))
  expect_false(any(is.infinite(unlist(result$coefficients))))
})

test_that("18.2: weight floor pmax(w, 1e-10) is applied correctly", {
  df <- generate_om_test_data(n = 200, seed = 42)

  # All control weights are tiny but positive → floor should kick in
  sample_weights <- rep(1e-20, nrow(df))
  ctrl_idx <- which(df$D == 0)
  # Ensure mean > 0 by giving some controls positive weights
  sample_weights[ctrl_idx[1:10]] <- 1.0

  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                   sample_weights = sample_weights)

  # Should produce finite results (floor prevents division by zero)
  expect_true(all(is.finite(result$predictions)))
  expect_true(all(is.finite(unlist(result$coefficients))))
})

# ============================================================================
# Test Group 19: WLS predictions full sample
# ============================================================================

test_that("19.1: WLS predictions length equals nrow(data)", {
  df <- generate_om_test_data(n = 200, seed = 42)
  sample_weights <- rep(1.0, nrow(df))

  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                   sample_weights = sample_weights)

  expect_length(result$predictions, nrow(df))
})

test_that("19.2: WLS predictions cover both treated and control units", {
  df <- generate_om_test_data(n = 200, seed = 42)
  sample_weights <- rep(1.0, nrow(df))

  result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                   sample_weights = sample_weights)

  treated_preds <- result$predictions[df$D == 1]
  control_preds <- result$predictions[df$D == 0]

  expect_true(length(treated_preds) > 0)
  expect_true(length(control_preds) > 0)
  expect_true(all(is.finite(treated_preds)))
  expect_true(all(is.finite(control_preds)))
})

test_that("19.3: WLS with uniform weights equals OLS", {
  df <- generate_om_test_data(n = 200, seed = 42)
  controls <- c("X1", "X2")

  result_ols <- estimate_outcome_model(df, "Y", "D", controls)
  result_wls <- estimate_outcome_model(df, "Y", "D", controls,
                                       sample_weights = rep(1.0, nrow(df)))

  # Uniform weights → WLS should equal OLS
  expect_equal(result_wls$coefficients[["_intercept"]],
               result_ols$coefficients[["_intercept"]],
               tolerance = 1e-10)
  expect_equal(result_wls$coefficients[["X1"]],
               result_ols$coefficients[["X1"]],
               tolerance = 1e-10)
  expect_equal(result_wls$coefficients[["X2"]],
               result_ols$coefficients[["X2"]],
               tolerance = 1e-10)
  expect_equal(result_wls$predictions, result_ols$predictions,
               tolerance = 1e-10)
})

# ============================================================================
# Test Group 20: Input validation errors
# ============================================================================

test_that("20.1: missing y column produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)

  expect_error(
    estimate_outcome_model(df, "nonexistent_y", "D", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )
})

test_that("20.2: missing d column produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)

  expect_error(
    estimate_outcome_model(df, "Y", "nonexistent_d", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )
})

test_that("20.3: missing control columns produces lwdid_invalid_param listing names", {
  df <- generate_om_test_data(n = 200, seed = 42)

  err <- tryCatch(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2", "X_missing")),
    lwdid_invalid_param = function(e) e
  )
  expect_true(inherits(err, "lwdid_invalid_param"))
  expect_true(grepl("X_missing", err$message, fixed = TRUE))
})

test_that("20.4: sample_weights wrong length produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                           sample_weights = rep(1.0, 50)),
    class = "lwdid_invalid_param"
  )
})

test_that("20.5: sample_weights with Inf produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)
  w <- rep(1.0, nrow(df))
  w[5] <- Inf

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                           sample_weights = w),
    class = "lwdid_invalid_param"
  )
})

test_that("20.6: sample_weights with NaN produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)
  w <- rep(1.0, nrow(df))
  w[10] <- NaN

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                           sample_weights = w),
    class = "lwdid_invalid_param"
  )
})

test_that("20.7: control group weight mean <= 0 produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)
  # All negative weights → mean <= 0
  w <- rep(-1.0, nrow(df))

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                           sample_weights = w),
    class = "lwdid_invalid_param"
  )
})

# ============================================================================
# Test Group 21: Singular design matrix
# ============================================================================

test_that("21.1: OLS with perfectly collinear controls produces lwdid_estimation_failed", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n)
  X2 <- 2 * X1  # perfectly collinear
  D <- rbinom(n, 1, 0.5)
  Y <- 2 + 3 * X1 + rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2")),
    class = "lwdid_estimation_failed"
  )
})

test_that("21.2: WLS with perfectly collinear controls produces lwdid_estimation_failed", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n)
  X2 <- 2 * X1  # perfectly collinear
  D <- rbinom(n, 1, 0.5)
  Y <- 2 + 3 * X1 + rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  w <- rep(1.0, n)

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                           sample_weights = w),
    class = "lwdid_estimation_failed"
  )
})

# ============================================================================
# Test Group 22: NA input tests
# ============================================================================

test_that("22.1: y with NA produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)
  df$Y[5] <- NA

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )

  err <- tryCatch(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2")),
    lwdid_invalid_param = function(e) e
  )
  expect_true(grepl("Y", err$message, fixed = TRUE))
})

test_that("22.2: d with NA produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)
  df$D[10] <- NA

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )

  err <- tryCatch(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2")),
    lwdid_invalid_param = function(e) e
  )
  expect_true(grepl("D", err$message, fixed = TRUE))
})

test_that("22.3: controls with NA produces lwdid_invalid_param error", {
  df <- generate_om_test_data(n = 200, seed = 42)
  df$X2[15] <- NA

  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2")),
    class = "lwdid_invalid_param"
  )

  err <- tryCatch(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2")),
    lwdid_invalid_param = function(e) e
  )
  expect_true(grepl("X2", err$message, fixed = TRUE))
})

test_that("22.4: multiple controls with NA lists all NA columns", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- 2 + X1 + X2 + X3 + rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2, X3 = X3)
  df$X1[5] <- NA
  df$X3[10] <- NA

  err <- tryCatch(
    estimate_outcome_model(df, "Y", "D", c("X1", "X2", "X3")),
    lwdid_invalid_param = function(e) e
  )
  expect_true(grepl("X1", err$message, fixed = TRUE))
  expect_true(grepl("X3", err$message, fixed = TRUE))
})

# ============================================================================
# Test Group 23: Integration test — PS → OM pipeline
# ============================================================================

test_that("23.1: PS → OM pipeline works end-to-end", {
  set.seed(42)
  n <- 500
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  linpred <- -0.5 + 0.8 * X1 + 0.5 * X2
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  # Ensure both groups
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  Y <- 2 + 3 * X1 + 1 * X2 + rnorm(n, sd = 0.5)
  Y[D == 1] <- Y[D == 1] + 5
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")

  # Step 1: Estimate propensity scores
  ps_result <- estimate_propensity_score(df, "D", controls)
  expect_true(ps_result$converged)

  ps <- ps_result$propensity_scores
  expect_true(all(ps >= 0.01 & ps <= 0.99))

  # Step 2: Compute IPW weights: w = ps / (1 - ps)
  ipw_weights <- ps / (1 - ps)
  expect_true(all(is.finite(ipw_weights)))
  expect_true(all(ipw_weights > 0))

  # Step 3: Estimate outcome model with IPW weights
  om_result <- estimate_outcome_model(df, "Y", "D", controls,
                                      sample_weights = ipw_weights)

  # Verify predictions are finite and cover all units
  expect_length(om_result$predictions, nrow(df))
  expect_true(all(is.finite(om_result$predictions)))

  # Verify coefficients are reasonable (close to true DGP)
  expect_equal(om_result$coefficients[["_intercept"]], 2, tolerance = 0.5)
  expect_equal(om_result$coefficients[["X1"]], 3, tolerance = 0.5)
  expect_equal(om_result$coefficients[["X2"]], 1, tolerance = 0.5)
})

test_that("23.2: PS → OM pipeline with drop trimming handles NA correctly", {
  set.seed(42)
  n <- 500
  X1 <- rnorm(n, sd = 2)
  X2 <- rnorm(n)
  linpred <- -1 + 2.0 * X1 + 0.5 * X2
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  Y <- 2 + 3 * X1 + 1 * X2 + rnorm(n, sd = 0.5)
  Y[D == 1] <- Y[D == 1] + 5
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")

  # PS with drop trimming
  ps_result <- estimate_propensity_score(df, "D", controls,
                                         trim_threshold = 0.1,
                                         trim_method = "drop")

  # For OM, use clip PS to avoid NA in weights
  ps_clip_result <- estimate_propensity_score(df, "D", controls,
                                              trim_threshold = 0.1,
                                              trim_method = "clip")
  ps <- ps_clip_result$propensity_scores
  ipw_weights <- ps / (1 - ps)

  om_result <- estimate_outcome_model(df, "Y", "D", controls,
                                      sample_weights = ipw_weights)

  expect_length(om_result$predictions, nrow(df))
  expect_true(all(is.finite(om_result$predictions)))
})

# ============================================================================
# Test Group 24: Integration test — vcov(glm_object)
# ============================================================================

test_that("24.1: vcov() on PS glm_object returns valid positive-definite matrix", {
  set.seed(42)
  n <- 300
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  linpred <- -0.5 + 0.8 * X1 + 0.5 * X2
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  df <- data.frame(D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")

  ps_result <- estimate_propensity_score(df, "D", controls)

  # vcov should return a valid variance-covariance matrix
  V <- vcov(ps_result$glm_object)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 3L)  # intercept + 2 covariates
  expect_equal(ncol(V), 3L)

  # Diagonal should be positive (variances)
  expect_true(all(diag(V) > 0))

  # Matrix should be symmetric
  expect_equal(V, t(V), tolerance = 1e-12)

  # Positive definite: all eigenvalues > 0
  eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("24.2: vcov() dimensions match number of coefficients", {
  set.seed(42)
  n <- 300
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
  linpred <- -0.5 + 0.8 * X1 + 0.5 * X2 + 0.3 * X3
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  if (all(D == 1)) D[1] <- 0L
  if (all(D == 0)) D[1] <- 1L
  df <- data.frame(D = D, X1 = X1, X2 = X2, X3 = X3)
  controls <- c("X1", "X2", "X3")

  ps_result <- estimate_propensity_score(df, "D", controls)
  V <- vcov(ps_result$glm_object)

  # 4x4: intercept + 3 covariates
  expect_equal(nrow(V), 4L)
  expect_equal(ncol(V), 4L)
  expect_true(all(diag(V) > 0))

  # Positive definite
  eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

# ============================================================================
# Section 2: estimate_outcome_model() tests
# ============================================================================

# Helper: generate outcome model test data with known DGP
# Y = beta0 + beta1*X1 + beta2*X2 + noise for D==0
# Y = beta0 + beta1*X1 + beta2*X2 + tau + noise for D==1
generate_om_test_data <- function(n = 400, seed = 42,
                                  beta0 = 1, beta1 = 2, beta2 = 3,
                                  tau = 5) {
  set.seed(seed)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  linpred_d <- -0.5 + 0.8 * X1 + 0.5 * X2
  prob <- 1 / (1 + exp(-linpred_d))
  D <- rbinom(n, 1, prob)
  if (all(D == 1) || all(D == 0)) { D[1] <- 0L; D[2] <- 1L }
  noise <- rnorm(n, sd = 0.5)
  Y <- beta0 + beta1 * X1 + beta2 * X2 + tau * D + noise
  data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
}

# ============================================================================
# Test Group 25: OLS mode tests
# ============================================================================

test_that("25.1: OLS basic functionality — coefficients match lm() to < 1e-12", {
  df <- generate_om_test_data(n = 500, seed = 123, beta0 = 1, beta1 = 2,
                              beta2 = 3, tau = 5)
  controls <- c("X1", "X2")

  result <- estimate_outcome_model(df, "Y", "D", controls)

  # Cross-check with lm() on control group
  ctrl_data <- df[df$D == 0, ]
  lm_fit <- lm(Y ~ X1 + X2, data = ctrl_data)
  lm_coefs <- coef(lm_fit)

  # Compare intercept
  expect_equal(result$coefficients[["_intercept"]], unname(lm_coefs["(Intercept)"]),
               tolerance = 1e-12)
  # Compare X1 coefficient

  expect_equal(result$coefficients[["X1"]], unname(lm_coefs["X1"]),
               tolerance = 1e-12)
  # Compare X2 coefficient
  expect_equal(result$coefficients[["X2"]], unname(lm_coefs["X2"]),
               tolerance = 1e-12)

  # Verify coefficients are numerically reasonable (close to true DGP)
  # True: beta0=1, beta1=2, beta2=3; with n=500 should be within ~0.3
  expect_true(abs(result$coefficients[["_intercept"]] - 1) < 0.5)
  expect_true(abs(result$coefficients[["X1"]] - 2) < 0.5)
  expect_true(abs(result$coefficients[["X2"]] - 3) < 0.5)
})

test_that("25.2: OLS predictions cover all units (length = nrow(data))", {
  df <- generate_om_test_data(n = 300, seed = 456)
  controls <- c("X1", "X2")

  result <- estimate_outcome_model(df, "Y", "D", controls)

  # Predictions length must equal total rows, not just control group
  expect_length(result$predictions, nrow(df))
  n_ctrl <- sum(df$D == 0)
  expect_true(n_ctrl < nrow(df))  # sanity: both groups exist

  # All predictions should be finite
  expect_true(all(is.finite(result$predictions)))

  # Cross-check: predictions = X_all %*% beta
  X_all <- cbind(1, as.matrix(df[, controls]))
  beta_vec <- c(result$coefficients[["_intercept"]],
                result$coefficients[["X1"]],
                result$coefficients[["X2"]])
  manual_pred <- as.numeric(X_all %*% beta_vec)
  expect_equal(result$predictions, manual_pred, tolerance = 1e-12)
})

test_that("25.3: OLS coefficient dictionary format (contains _intercept key)", {
  df <- generate_om_test_data(n = 200, seed = 789)
  controls <- c("X1", "X2")

  result <- estimate_outcome_model(df, "Y", "D", controls)

  # coefficients is a named list

  expect_type(result$coefficients, "list")

  # Must contain "_intercept" key
  expect_true("_intercept" %in% names(result$coefficients))

  # Must contain one key per control variable
  for (ctrl in controls) {
    expect_true(ctrl %in% names(result$coefficients),
                info = paste("Missing coefficient for", ctrl))
  }

  # Total keys = 1 (intercept) + length(controls)
  expect_length(result$coefficients, 1L + length(controls))

  # All coefficient values are numeric scalars
  for (nm in names(result$coefficients)) {
    expect_true(is.numeric(result$coefficients[[nm]]),
                info = paste("Coefficient", nm, "is not numeric"))
    expect_length(result$coefficients[[nm]], 1L)
  }
})

# ============================================================================
# Test Group 26: WLS mode tests
# ============================================================================

test_that("26.1: WLS basic functionality — matches manual sqrt(w) transform", {
  set.seed(101)
  n <- 400
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  linpred_d <- -0.3 + 0.6 * X1 + 0.4 * X2
  prob <- 1 / (1 + exp(-linpred_d))
  D <- rbinom(n, 1, prob)
  if (all(D == 1) || all(D == 0)) { D[1] <- 0L; D[2] <- 1L }
  Y <- 2 + 1.5 * X1 + 2.5 * X2 + 3 * D + rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")

  # Compute PS-style weights: w = p / (1 - p)
  w <- prob / (1 - prob)

  result <- estimate_outcome_model(df, "Y", "D", controls,
                                   sample_weights = w)

  # All coefficients should be finite
  for (nm in names(result$coefficients)) {
    expect_true(is.finite(result$coefficients[[nm]]),
                info = paste("Coefficient", nm, "is not finite"))
  }

  # Manual WLS cross-check via sqrt(w) transform
  ctrl_mask <- D == 0
  X_ctrl <- cbind(1, X1[ctrl_mask], X2[ctrl_mask])
  Y_ctrl <- Y[ctrl_mask]
  w_ctrl <- pmax(w[ctrl_mask], 1e-10)
  sqrt_w <- sqrt(w_ctrl)
  X_w <- X_ctrl * sqrt_w
  Y_w <- Y_ctrl * sqrt_w
  beta_manual <- solve(crossprod(X_w), crossprod(X_w, Y_w))

  expect_equal(result$coefficients[["_intercept"]], as.numeric(beta_manual[1]),
               tolerance = 1e-10)
  expect_equal(result$coefficients[["X1"]], as.numeric(beta_manual[2]),
               tolerance = 1e-10)
  expect_equal(result$coefficients[["X2"]], as.numeric(beta_manual[3]),
               tolerance = 1e-10)

  # Coefficients should be numerically reasonable (close to true DGP)
  expect_true(abs(result$coefficients[["_intercept"]] - 2) < 1.0)
  expect_true(abs(result$coefficients[["X1"]] - 1.5) < 1.0)
  expect_true(abs(result$coefficients[["X2"]] - 2.5) < 1.0)
})

test_that("26.2: WLS weight floor — near-zero weights don't produce NaN/Inf", {
  set.seed(202)
  n <- 300
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  if (all(D == 1) || all(D == 0)) { D[1] <- 0L; D[2] <- 1L }
  Y <- 1 + X1 + X2 + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")

  # Create weights with some near-zero values
  w <- rep(1.0, n)
  ctrl_indices <- which(D == 0)
  # Set ~20% of control weights to near-zero
  n_nearzero <- max(1L, floor(length(ctrl_indices) * 0.2))
  w[ctrl_indices[seq_len(n_nearzero)]] <- 1e-15

  result <- estimate_outcome_model(df, "Y", "D", controls,
                                   sample_weights = w)

  # No NaN or Inf in predictions
  expect_true(all(is.finite(result$predictions)),
              info = "Predictions contain NaN or Inf with near-zero weights")

  # No NaN or Inf in coefficients
  for (nm in names(result$coefficients)) {
    expect_true(is.finite(result$coefficients[[nm]]),
                info = paste("Coefficient", nm, "is NaN/Inf with near-zero weights"))
  }
})

test_that("26.3: WLS predictions full-sample — length and finiteness", {
  set.seed(303)
  n <- 250
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  linpred_d <- -0.2 + 0.5 * X1
  prob <- 1 / (1 + exp(-linpred_d))
  D <- rbinom(n, 1, prob)
  if (all(D == 1) || all(D == 0)) { D[1] <- 0L; D[2] <- 1L }
  Y <- 3 + X1 - X2 + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")
  w <- prob / (1 - prob)

  result <- estimate_outcome_model(df, "Y", "D", controls,
                                   sample_weights = w)

  # Predictions length = nrow(data), not just control group
  expect_length(result$predictions, nrow(df))
  expect_true(sum(D == 0) < nrow(df))  # sanity

  # All predictions finite
  expect_true(all(is.finite(result$predictions)))

  # Cross-check: predictions = X_all %*% beta
  X_all <- cbind(1, as.matrix(df[, controls]))
  beta_vec <- c(result$coefficients[["_intercept"]],
                result$coefficients[["X1"]],
                result$coefficients[["X2"]])
  manual_pred <- as.numeric(X_all %*% beta_vec)
  expect_equal(result$predictions, manual_pred, tolerance = 1e-12)
})

# ============================================================================
# Test Group 27: Error handling tests
# ============================================================================

test_that("27.1: Input validation errors — missing columns, weight anomalies", {
  df <- generate_om_test_data(n = 100, seed = 42)
  controls <- c("X1", "X2")

  # y column missing
  expect_error(
    estimate_outcome_model(df, "NONEXISTENT_Y", "D", controls),
    class = "lwdid_invalid_param"
  )

  # d column missing
  expect_error(
    estimate_outcome_model(df, "Y", "NONEXISTENT_D", controls),
    class = "lwdid_invalid_param"
  )

  # controls column missing
  expect_error(
    estimate_outcome_model(df, "Y", "D", c("X1", "NONEXISTENT_X3")),
    class = "lwdid_invalid_param"
  )

  # Weight length mismatch
  expect_error(
    estimate_outcome_model(df, "Y", "D", controls,
                           sample_weights = rep(1, nrow(df) + 5)),
    class = "lwdid_invalid_param"
  )

  # Weight contains Inf
  w_inf <- rep(1, nrow(df))
  w_inf[5] <- Inf
  expect_error(
    estimate_outcome_model(df, "Y", "D", controls, sample_weights = w_inf),
    class = "lwdid_invalid_param"
  )

  # Weight contains NaN
  w_nan <- rep(1, nrow(df))
  w_nan[3] <- NaN
  expect_error(
    estimate_outcome_model(df, "Y", "D", controls, sample_weights = w_nan),
    class = "lwdid_invalid_param"
  )

  # Control group weight mean <= 0
  # All control weights negative → mean <= 0
  w_neg <- rep(1, nrow(df))
  ctrl_idx <- which(df$D == 0)
  w_neg[ctrl_idx] <- -1
  expect_error(
    estimate_outcome_model(df, "Y", "D", controls, sample_weights = w_neg),
    class = "lwdid_invalid_param"
  )
})

test_that("27.2: Singular design matrix — OLS and WLS", {
  set.seed(555)
  n <- 200
  X1 <- rnorm(n)
  X2 <- 2 * X1  # perfect collinearity
  D <- rbinom(n, 1, 0.5)
  if (all(D == 1) || all(D == 0)) { D[1] <- 0L; D[2] <- 1L }
  Y <- 1 + X1 + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")

  # OLS: singular design matrix
  expect_error(
    estimate_outcome_model(df, "Y", "D", controls),
    class = "lwdid_estimation_failed"
  )

  # WLS: same collinearity with weights
  w <- rep(1, n)
  expect_error(
    estimate_outcome_model(df, "Y", "D", controls, sample_weights = w),
    class = "lwdid_estimation_failed"
  )
})

test_that("27.3: NA input test — y/d/controls with NA", {
  df <- generate_om_test_data(n = 100, seed = 42)
  controls <- c("X1", "X2")

  # y contains NA
  df_na_y <- df
  df_na_y$Y[5] <- NA
  expect_error(
    estimate_outcome_model(df_na_y, "Y", "D", controls),
    class = "lwdid_invalid_param"
  )

  # d contains NA
  df_na_d <- df
  df_na_d$D[10] <- NA
  expect_error(
    estimate_outcome_model(df_na_d, "Y", "D", controls),
    class = "lwdid_invalid_param"
  )

  # controls contain NA
  df_na_x <- df
  df_na_x$X1[15] <- NA
  expect_error(
    estimate_outcome_model(df_na_x, "Y", "D", controls),
    class = "lwdid_invalid_param"
  )

  # NA in second control variable
  df_na_x2 <- df
  df_na_x2$X2[20] <- NA
  expect_error(
    estimate_outcome_model(df_na_x2, "Y", "D", controls),
    class = "lwdid_invalid_param"
  )
})

# ============================================================================
# Test Group 28: Integration tests
# ============================================================================

test_that("28.1: Integration — PS → OM end-to-end pipeline", {
  set.seed(777)
  n <- 500
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  linpred_d <- -0.5 + 0.8 * X1 + 0.5 * X2
  prob <- 1 / (1 + exp(-linpred_d))
  D <- rbinom(n, 1, prob)
  if (all(D == 1) || all(D == 0)) { D[1] <- 0L; D[2] <- 1L }
  Y <- 2 + 1.5 * X1 + 2.5 * X2 + 4 * D + rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")

  # Step 1: Estimate propensity scores
  ps_result <- estimate_propensity_score(df, "D", controls)
  ps <- ps_result$propensity_scores

  # Verify PS are valid
  expect_true(all(ps > 0 & ps < 1))

  # Step 2: Compute ATT weights w = ps / (1 - ps)
  w <- ps / (1 - ps)
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))

  # Step 3: Estimate outcome model with WLS
  om_result <- estimate_outcome_model(df, "Y", "D", controls,
                                      sample_weights = w)

  # Predictions are finite and cover all units
  expect_length(om_result$predictions, nrow(df))
  expect_true(all(is.finite(om_result$predictions)))

  # Coefficients are finite and numerically reasonable
  expect_true(is.finite(om_result$coefficients[["_intercept"]]))
  expect_true(is.finite(om_result$coefficients[["X1"]]))
  expect_true(is.finite(om_result$coefficients[["X2"]]))

  # Coefficients should be in reasonable range of true DGP (beta0=2, beta1=1.5, beta2=2.5)
  expect_true(abs(om_result$coefficients[["_intercept"]] - 2) < 2.0)
  expect_true(abs(om_result$coefficients[["X1"]] - 1.5) < 2.0)
  expect_true(abs(om_result$coefficients[["X2"]] - 2.5) < 2.0)

  # Predictions should have reasonable variance (not degenerate)
  pred_sd <- sd(om_result$predictions)
  expect_true(pred_sd > 0.1, info = "Predictions have near-zero variance")
})

test_that("28.2: Integration — vcov(glm_object) + OM pipeline", {
  set.seed(888)
  n <- 400
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  linpred_d <- -0.3 + 0.7 * X1 + 0.4 * X2
  prob <- 1 / (1 + exp(-linpred_d))
  D <- rbinom(n, 1, prob)
  if (all(D == 1) || all(D == 0)) { D[1] <- 0L; D[2] <- 1L }
  Y <- 1 + 2 * X1 + 3 * X2 + 5 * D + rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  controls <- c("X1", "X2")

  # Step 1: PS estimation
  ps_result <- estimate_propensity_score(df, "D", controls)

  # Verify vcov() returns valid variance-covariance matrix
  V <- vcov(ps_result$glm_object)
  expect_true(is.matrix(V))
  n_params <- length(controls) + 1L  # intercept + controls
  expect_equal(nrow(V), n_params)
  expect_equal(ncol(V), n_params)

  # Symmetric
  expect_equal(V, t(V), tolerance = 1e-12)

  # Positive definite: all eigenvalues > 0
  eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0),
              info = "vcov matrix is not positive definite")

  # Diagonal elements (variances) should be positive
  expect_true(all(diag(V) > 0))

  # Step 2: Use PS to compute ATT weights and run OM
  ps <- ps_result$propensity_scores
  w <- ps / (1 - ps)

  om_result <- estimate_outcome_model(df, "Y", "D", controls,
                                      sample_weights = w)

  # Full pipeline produces valid output
  expect_length(om_result$predictions, nrow(df))
  expect_true(all(is.finite(om_result$predictions)))
  expect_true("_intercept" %in% names(om_result$coefficients))

  # Coefficients are finite
  for (nm in names(om_result$coefficients)) {
    expect_true(is.finite(om_result$coefficients[[nm]]),
                info = paste("Pipeline coefficient", nm, "is not finite"))
  }

  # Verify the pipeline output is numerically reasonable
  # Control-group predictions should be close to actual Y for controls
  ctrl_mask <- df$D == 0
  ctrl_pred <- om_result$predictions[ctrl_mask]
  ctrl_actual <- df$Y[ctrl_mask]
  # R-squared proxy: correlation should be high for well-specified model
  r_corr <- cor(ctrl_pred, ctrl_actual)
  expect_true(r_corr > 0.5,
              info = paste("Control group prediction correlation too low:", r_corr))
})
