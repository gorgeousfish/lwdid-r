library(testthat)

if (!exists("lwdid", mode = "function")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to run test-monte-carlo-e606.R directly.")
  }

  devtools::load_all(
    "/Users/cxy/Desktop/lwdid_r/lwdid-r",
    export_all = FALSE,
    quiet = TRUE
  )
}

quiet_e606_monte_carlo <- function(expr) {
  withCallingHandlers(
    expr,
    lwdid_small_sample = function(w) invokeRestart("muffleWarning"),
    lwdid_overlap = function(w) invokeRestart("muffleWarning"),
    lwdid_data = function(w) invokeRestart("muffleWarning"),
    lwdid_numerical = function(w) invokeRestart("muffleWarning"),
    lwdid_convergence = function(w) invokeRestart("muffleWarning"),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

e606_monte_carlo_reps <- function() {
  reps <- as.integer(Sys.getenv("LWDID_E606_MONTE_CARLO_REPS", "1000"))
  if (is.na(reps) || reps < 50L) {
    stop("LWDID_E606_MONTE_CARLO_REPS must be an integer >= 50.")
  }
  reps
}

e606_monte_carlo_cache <- new.env(parent = emptyenv())

build_e606_monte_carlo_sample <- function(seed) {
  set.seed(seed)

  n <- 200L
  x <- stats::rnorm(n)
  treat_prob <- stats::plogis(0.5 * x)
  d <- stats::rbinom(n, 1L, treat_prob)

  while (sum(d) < 20L || sum(d == 0L) < 20L) {
    d <- stats::rbinom(n, 1L, treat_prob)
  }

  y0 <- x + stats::rnorm(n)
  y1 <- y0 + 1.0

  data.frame(
    Y = ifelse(d == 1L, y1, y0),
    D = d,
    X = x
  )
}

build_e606_ipwra_dr_sample <- function(seed) {
  sample_df <- build_e606_monte_carlo_sample(seed)
  sample_df$X_sq <- sample_df$X^2
  sample_df
}

run_e606_monte_carlo <- function(estimator,
                                 reps = e606_monte_carlo_reps(),
                                 base_seed = 42L,
                                 psm_se_method = "abadie_imbens") {
  cache_key <- paste(estimator, reps, base_seed, psm_se_method, sep = "::")
  if (exists(cache_key, envir = e606_monte_carlo_cache, inherits = FALSE)) {
    return(get(cache_key, envir = e606_monte_carlo_cache, inherits = FALSE))
  }

  att <- numeric(reps)
  se <- numeric(reps)
  covered <- logical(reps)

  for (b in seq_len(reps)) {
    sample_df <- build_e606_monte_carlo_sample(base_seed + b - 1L)

    fit <- quiet_e606_monte_carlo(
      if (identical(estimator, "ipw")) {
        lwdid:::estimate_ipw(sample_df, "Y", "D", "X")
      } else if (identical(estimator, "psm")) {
        lwdid::estimate_psm(
          sample_df,
          "Y",
          "D",
          "X",
          n_neighbors = 1L,
          with_replacement = TRUE,
          se_method = psm_se_method
        )
      } else {
        lwdid:::estimate_ipwra(
          sample_df,
          "Y",
          "D",
          controls = "X",
          propensity_controls = "X"
        )
      }
    )

    att[[b]] <- fit$att
    se[[b]] <- fit$se
    covered[[b]] <- fit$ci_lower <= 1.0 && 1.0 <= fit$ci_upper
  }

  result <- list(
    reps = reps,
    mean_att = mean(att),
    mean_bias = mean(att - 1.0),
    mean_se = mean(se),
    coverage = mean(covered),
    all_finite = all(is.finite(att)) && all(is.finite(se))
  )

  assign(cache_key, result, envir = e606_monte_carlo_cache)
  result
}

run_e606_ipwra_dr_monte_carlo <- function(
  scenario,
  reps = e606_monte_carlo_reps(),
  base_seed = 42L
) {
  if (!scenario %in% c("ps_misspecified", "outcome_misspecified")) {
    stop("scenario must be one of 'ps_misspecified' or 'outcome_misspecified'.")
  }

  controls <- if (identical(scenario, "ps_misspecified")) "X" else "X_sq"
  propensity_controls <- if (identical(scenario, "ps_misspecified")) "X_sq" else "X"

  att <- numeric(reps)
  se <- numeric(reps)
  covered <- logical(reps)

  for (b in seq_len(reps)) {
    sample_df <- build_e606_ipwra_dr_sample(base_seed + b - 1L)

    fit <- quiet_e606_monte_carlo(
      lwdid:::estimate_ipwra(
        sample_df,
        "Y",
        "D",
        controls = controls,
        propensity_controls = propensity_controls
      )
    )

    att[[b]] <- fit$att
    se[[b]] <- fit$se
    covered[[b]] <- fit$ci_lower <= 1.0 && 1.0 <= fit$ci_upper
  }

  list(
    scenario = scenario,
    reps = reps,
    mean_att = mean(att),
    mean_bias = mean(att - 1.0),
    mean_se = mean(se),
    coverage = mean(covered),
    all_finite = all(is.finite(att)) && all(is.finite(se))
  )
}

test_that("TC-6.6.33: IPW Monte Carlo bias stays below the acceptance band under correct PS specification", {
  skip_on_cran()

  summary <- run_e606_monte_carlo("ipw")

  expect_equal(summary$reps, e606_monte_carlo_reps())
  expect_true(summary$all_finite)
  expect_gt(summary$coverage, 0.85)
  expect_lte(abs(summary$mean_bias), 0.1)
})

test_that("TC-6.6.34: IPWRA Monte Carlo bias stays below the tighter doubly-robust band", {
  skip_on_cran()

  summary <- run_e606_monte_carlo("ipwra")

  expect_equal(summary$reps, e606_monte_carlo_reps())
  expect_true(summary$all_finite)
  expect_gt(summary$mean_se, 0)
  expect_lte(abs(summary$mean_bias), 0.05)
})

test_that("TC-6.6.35: IPWRA analytical coverage stays within the acceptance band", {
  skip_on_cran()

  summary <- run_e606_monte_carlo("ipwra")

  expect_equal(summary$reps, e606_monte_carlo_reps())
  expect_true(summary$all_finite)
  expect_gte(summary$coverage, 0.93)
  expect_lte(summary$coverage, 0.97)
})

test_that("TC-6.6.36: IPWRA stays more efficient than IPW on the shared DGP", {
  skip_on_cran()

  ipw_summary <- run_e606_monte_carlo("ipw")
  ipwra_summary <- run_e606_monte_carlo("ipwra")

  expect_equal(ipw_summary$reps, e606_monte_carlo_reps())
  expect_equal(ipwra_summary$reps, e606_monte_carlo_reps())
  expect_true(ipw_summary$all_finite)
  expect_true(ipwra_summary$all_finite)
  expect_lt(ipwra_summary$mean_se, ipw_summary$mean_se)
})

test_that("TC-6.6.37: IPWRA stays near the true ATT when the propensity model is misspecified but the outcome model is correct", {
  skip_on_cran()

  summary <- run_e606_ipwra_dr_monte_carlo("ps_misspecified")

  expect_equal(summary$reps, e606_monte_carlo_reps())
  expect_identical(summary$scenario, "ps_misspecified")
  expect_true(summary$all_finite)
  expect_lte(abs(summary$mean_bias), 0.3)
})

test_that("TC-6.6.38: IPWRA stays near the true ATT when the propensity model is correct but the outcome model is misspecified", {
  skip_on_cran()

  summary <- run_e606_ipwra_dr_monte_carlo("outcome_misspecified")

  expect_equal(summary$reps, e606_monte_carlo_reps())
  expect_identical(summary$scenario, "outcome_misspecified")
  expect_true(summary$all_finite)
  expect_lte(abs(summary$mean_bias), 0.3)
})

test_that("TC-6.6.39: PSM bias stays within the acceptance band on the shared DGP", {
  skip_on_cran()

  summary <- run_e606_monte_carlo("psm")

  expect_equal(summary$reps, e606_monte_carlo_reps())
  expect_true(summary$all_finite)
  expect_lte(abs(summary$mean_bias), 0.15)
})

test_that("TC-6.6.40: PSM coverage stays within the acceptance band under the reuse-adjusted CI path", {
  skip_on_cran()

  summary <- run_e606_monte_carlo("psm", psm_se_method = "reuse_adjusted")

  expect_equal(summary$reps, e606_monte_carlo_reps())
  expect_true(summary$all_finite)
  expect_gte(summary$coverage, 0.92)
  expect_lte(summary$coverage, 0.99)
})

test_that("TC-6.6.41: IPW analytical coverage stays within the acceptance band", {
  skip_on_cran()

  summary <- run_e606_monte_carlo("ipw")

  expect_equal(summary$reps, e606_monte_carlo_reps())
  expect_true(summary$all_finite)
  expect_gte(summary$coverage, 0.93)
  expect_lte(summary$coverage, 0.97)
})
