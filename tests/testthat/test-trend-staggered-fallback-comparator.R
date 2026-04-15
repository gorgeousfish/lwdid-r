library(testthat)

resolve_staggered_fallback_oracle_path <- function() {
  file.path(
    "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity",
    "20260325-qa-parity-e8-06-trend-public-staggered-fallback.json"
  )
}

read_staggered_fallback_oracle <- function() {
  jsonlite::fromJSON(
    resolve_staggered_fallback_oracle_path(),
    simplifyVector = FALSE
  )
}

test_that("E8-06 staggered fallback oracle normalizes no-estimate joint F as null", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_staggered_fallback_oracle_path()
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend staggered fallback oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  case <- read_staggered_fallback_oracle()$cases$staggered_simple_did_fallback_public_api

  expect_null(case$python_result$joint_f_stat)
  expect_null(case$r_current_behavior$joint_f_stat)
  expect_null(case$comparison$python_joint_f_stat)
  expect_null(case$comparison$r_joint_f_stat)
})
