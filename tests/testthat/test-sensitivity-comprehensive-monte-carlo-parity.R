library(testthat)

resolve_monte_carlo_parity_path <- function(filename) {
  resolve_parity_fixture_path(filename)
}

compute_current_monte_carlo_summary <- function(fixture_df) {
  split_fixture <- split(fixture_df, fixture_df$replication)

  rows <- lapply(split_fixture, function(rep_df) {
    true_att <- unique(rep_df$true_att)
    expect_length(true_att, 1L)

    warnings_seen <- character(0)
    result <- suppressMessages(
      withCallingHandlers(
        lwdid_sensitivity(
          data = rep_df[c("id", "year", "y", "d", "post", "x1", "x2")],
          y = "y",
          ivar = "id",
          tvar = "year",
          d = "d",
          post = "post",
          controls = c("x1", "x2"),
          type = "all",
          verbose = FALSE
        ),
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    )

    demean_att <- result$transformation_comparison$demean_att
    demean_se <- result$transformation_comparison$demean_se
    ci_lower <- demean_att - 1.96 * demean_se
    ci_upper <- demean_att + 1.96 * demean_se

    data.frame(
      replication = unique(rep_df$replication),
      true_att = true_att,
      demean_att = demean_att,
      demean_se = demean_se,
      bias = demean_att - true_att,
      covered = ci_lower <= true_att && true_att <= ci_upper
    )
  })

  metrics <- do.call(rbind, rows)
  metrics$covered <- as.logical(metrics$covered)

  list(
    replications = nrow(metrics),
    summary = list(
      mean_true_att = mean(metrics$true_att),
      mean_att = mean(metrics$demean_att),
      mean_se = mean(metrics$demean_se),
      mean_bias = mean(metrics$bias),
      mean_abs_bias = mean(abs(metrics$bias)),
      coverage = mean(metrics$covered),
      finite_att = all(is.finite(metrics$demean_att)),
      finite_se = all(is.finite(metrics$demean_se))
    ),
    by_replication = metrics
  )
}

test_that("E8-03 Monte Carlo oracle matches current common-timing sensitivity metrics", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_monte_carlo_parity_path(
    "20260323-qa-parity-e8-03-monte-carlo-common-timing.json"
  )
  fixture_path <- resolve_monte_carlo_parity_path(
    "e8_03_monte_carlo_common_timing_fixture.csv"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing Monte Carlo parity oracle:", oracle_path)
  )
  expect_true(
    file.exists(fixture_path),
    info = paste("missing Monte Carlo fixture:", fixture_path)
  )
  if (!file.exists(oracle_path) || !file.exists(fixture_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  fixture_df <- read.csv(fixture_path, stringsAsFactors = FALSE)
  current <- compute_current_monte_carlo_summary(fixture_df)

  expect_true(isTRUE(oracle$comparison$status == "passed"))
  expect_true(isTRUE(oracle$comparison$numeric_status == "passed"))
  expect_true(current$summary$finite_att)
  expect_true(current$summary$finite_se)

  py_ref <- oracle$python_reference$summary

  expect_equal(
    current$summary$mean_att,
    py_ref$mean_att,
    tolerance = 1e-6
  )
  expect_equal(
    current$summary$mean_se,
    py_ref$mean_se,
    tolerance = 1e-4
  )
  expect_equal(
    current$summary$mean_bias,
    py_ref$mean_bias,
    tolerance = 1e-6
  )
  expect_equal(
    current$summary$mean_abs_bias,
    py_ref$mean_abs_bias,
    tolerance = 1e-6
  )
  expect_equal(
    current$summary$coverage,
    py_ref$coverage,
    tolerance = 1e-12
  )

  tolerances <- oracle$tolerances
  expect_lte(abs(current$summary$mean_bias), tolerances$monte_carlo_bias)
  expect_lte(abs(current$summary$coverage - 0.95), tolerances$monte_carlo_coverage)
})
