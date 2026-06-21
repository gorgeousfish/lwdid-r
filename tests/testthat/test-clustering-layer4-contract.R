test_that("E8-04 Layer 4 clustering contract parses completely", {
  skip_if_not_installed("yaml")

  contract_path <- test_path(
    "_fixtures",
    "parity",
    "e8_04_clustering_layer4_monte_carlo_contract.yaml"
  )
  expect_true(file.exists(contract_path))

  contract <- yaml::yaml.load(
    paste(readLines(contract_path, encoding = "UTF-8", warn = FALSE),
          collapse = "\n")
  )

  expected_ids <- c(
    "coverage_large_g",
    "coverage_medium_g",
    "coverage_small_g",
    "wild_bootstrap_coverage_small_g",
    "wild_bootstrap_relative_gap",
    "size_under_null",
    "power_under_alternative",
    "unbalanced_treatment_coverage",
    "heterogeneous_cluster_sizes_coverage"
  )
  scenario_ids <- vapply(
    contract$scenarios,
    function(scenario) scenario$id,
    character(1)
  )

  expect_identical(scenario_ids, expected_ids)
  expect_identical(anyDuplicated(scenario_ids), 0L)
  expect_identical(
    contract$scenarios[[4L]]$actual_n_bootstrap,
    1024L
  )
  expect_identical(
    contract$scenarios[[5L]]$estimators$wild$actual_n_bootstrap,
    1024L
  )
})
