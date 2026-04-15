library(testthat)

resolve_clustering_parity_path <- function(filename) {
  resolved <- resolve_parity_fixture_path(filename)
  if (file.exists(resolved)) {
    return(resolved)
  }

  testthat::skip(
    paste(
      "clustering parity artifacts are unavailable in this test environment:",
      resolved
    )
  )
}

read_clustering_fixture <- function(filename) {
  utils::read.csv(
    resolve_clustering_parity_path(filename),
    stringsAsFactors = FALSE
  )
}

read_clustering_oracle <- function(filename) {
  jsonlite::fromJSON(
    resolve_clustering_parity_path(filename),
    simplifyVector = FALSE
  )
}

as_character_vector <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }

  vec <- unlist(x, use.names = FALSE)
  if (length(vec) == 0L) {
    return(character(0))
  }

  as.character(vec)
}

as_named_int_vector <- function(x) {
  if (is.null(x)) {
    return(integer(0))
  }

  vec <- as.integer(unlist(x, use.names = TRUE))
  stats::setNames(vec, names(unlist(x, use.names = TRUE)))
}

expect_cluster_stats_match <- function(actual, expected, tolerance = 1e-10) {
  expected_cluster_sizes <- as_named_int_vector(expected$cluster_sizes)
  actual_cluster_sizes <- actual$cluster_sizes

  actual_cluster_sizes <- actual_cluster_sizes[order(names(actual_cluster_sizes))]
  expected_cluster_sizes <- expected_cluster_sizes[order(names(expected_cluster_sizes))]

  expect_identical(actual$var_name, expected$var_name)
  expect_identical(actual$n_clusters, expected$n_clusters)
  expect_equal(actual_cluster_sizes, expected_cluster_sizes)
  expect_identical(actual$min_size, expected$min_size)
  expect_identical(actual$max_size, expected$max_size)
  expect_equal(actual$mean_size, expected$mean_size, tolerance = tolerance)
  expect_equal(actual$median_size, expected$median_size, tolerance = tolerance)
  expect_equal(actual$cv, expected$cv, tolerance = tolerance)
  expect_identical(actual$n_treated_clusters, expected$n_treated_clusters)
  expect_identical(actual$n_control_clusters, expected$n_control_clusters)
  expect_identical(actual$treatment_varies_within, expected$treatment_varies_within)
  expect_identical(
    actual$n_clusters_with_variation,
    expected$n_clusters_with_variation
  )
  expect_identical(actual$is_nested_in_unit, expected$is_nested_in_unit)
  expect_identical(actual$level_relative_to_unit, expected$level_relative_to_unit)
  expect_equal(actual$units_per_cluster, expected$units_per_cluster, tolerance = tolerance)
  expect_equal(actual$balance_score, expected$balance_score, tolerance = tolerance)
  expect_equal(actual$cv_score, expected$cv_score, tolerance = tolerance)
  expect_equal(
    actual$reliability_score,
    expected$reliability_score,
    tolerance = tolerance
  )
  expect_identical(actual$is_valid, expected$is_valid)
}

compute_cluster_unique_counts <- function(data, cluster_var, treat_var) {
  counts <- data.table::as.data.table(data)[
    ,
    .(n_unique = data.table::uniqueN(get(treat_var), na.rm = TRUE)),
    by = cluster_var
  ]

  stats::setNames(
    as.integer(counts$n_unique),
    as.character(counts[[cluster_var]])
  )
}

test_that("E8-04 Task 1 parity oracle exists and stays aligned on hierarchical fixture", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_clustering_parity_path(
    "20260324-qa-parity-e8-04-clustering-task1-python-oracle.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing clustering parity oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  case <- oracle$cases$hierarchical

  fixture <- read_clustering_fixture(case$fixture_csv)
  expected_state <- case$analyze_cluster_var$state
  expected_county <- case$analyze_cluster_var$county

  actual_state <- .analyze_cluster_var(
    fixture,
    ivar = case$ivar,
    cluster_var = "state",
    gvar = case$gvar,
    d = NULL
  )
  actual_county <- .analyze_cluster_var(
    fixture,
    ivar = case$ivar,
    cluster_var = "county",
    gvar = case$gvar,
    d = NULL
  )

  expect_identical(
    .determine_cluster_level(fixture, case$ivar, "state"),
    case$cluster_levels$state
  )
  expect_identical(
    .determine_cluster_level(fixture, case$ivar, "county"),
    case$cluster_levels$county
  )
  expect_cluster_stats_match(actual_state, expected_state)
  expect_cluster_stats_match(actual_county, expected_county)
})

test_that("E8-04 Task 1 parity oracle covers within-cluster variation and never-treated masks", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_clustering_parity_path(
    "20260324-qa-parity-e8-04-clustering-task1-python-oracle.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing clustering parity oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)

  variation_case <- oracle$cases$within_cluster_variation
  variation_fixture <- read_clustering_fixture(variation_case$fixture_csv)
  actual_variation <- .analyze_cluster_var(
    variation_fixture,
    ivar = variation_case$ivar,
    cluster_var = "state",
    gvar = variation_case$gvar,
    d = NULL
  )
  expect_cluster_stats_match(
    actual_variation,
    variation_case$analyze_cluster_var$state
  )

  never_treated_case <- oracle$cases$never_treated
  never_treated_fixture <- read_clustering_fixture(never_treated_case$fixture_csv)
  actual_never_treated <- .analyze_cluster_var(
    never_treated_fixture,
    ivar = never_treated_case$ivar,
    cluster_var = "cluster",
    gvar = never_treated_case$gvar,
    d = NULL
  )
  expect_cluster_stats_match(
    actual_never_treated,
    never_treated_case$analyze_cluster_var$cluster
  )
})

test_that("E8-04 Task 3 public-contract oracle freezes small-cluster diagnosis semantics", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_clustering_parity_path(
    "20260324-qa-parity-e8-04-clustering-task23-contract.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing clustering public-contract oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_clustering_oracle(
    "20260324-qa-parity-e8-04-clustering-task23-contract.json"
  )
  case <- oracle$cases$small_cluster
  fixture <- read_clustering_fixture(case$fixture_csv)

  actual_diag <- diagnose_clustering(
    fixture,
    ivar = case$ivar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = case$gvar,
    d = NULL,
    verbose = FALSE
  )

  expect_identical(
    actual_diag$recommended_cluster_var,
    case$python_diagnose_clustering$recommended_cluster_var
  )
  expect_identical(
    actual_diag$treatment_variation_level,
    case$python_diagnose_clustering$treatment_variation_level
  )
  expect_identical(
    actual_diag$recommendation_reason,
    case$python_diagnose_clustering$recommendation_reason
  )
  expect_identical(
    actual_diag$warnings,
    as_character_vector(case$python_diagnose_clustering$warnings)
  )

  expect_identical(
    case$python_recommend_clustering_level$recommended_var,
    "state"
  )
  expect_identical(
    case$python_recommend_clustering_level$n_clusters,
    5L
  )
  expect_true(
    isTRUE(case$python_recommend_clustering_level$use_wild_bootstrap)
  )

  actual_rec <- recommend_clustering(
    fixture,
    ivar = case$ivar,
    tvar = "year",
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = case$gvar,
    d = NULL,
    verbose = FALSE
  )

  expect_identical(
    actual_rec$recommended_var,
    case$python_recommend_clustering_level$recommended_var
  )
  expect_identical(
    actual_rec$n_clusters,
    case$python_recommend_clustering_level$n_clusters
  )
  expect_identical(
    actual_rec$n_treated_clusters,
    case$python_recommend_clustering_level$n_treated_clusters
  )
  expect_identical(
    actual_rec$n_control_clusters,
    case$python_recommend_clustering_level$n_control_clusters
  )
  expect_equal(
    actual_rec$confidence,
    case$python_recommend_clustering_level$confidence,
    tolerance = 1e-12
  )
  expect_identical(
    actual_rec$reasons,
    unlist(case$python_recommend_clustering_level$reasons, use.names = FALSE)
  )
  expect_identical(
    actual_rec$warnings,
    unlist(case$python_recommend_clustering_level$warnings, use.names = FALSE)
  )
  expect_identical(
    actual_rec$wild_bootstrap_reason,
    case$python_recommend_clustering_level$wild_bootstrap_reason
  )
  expect_true(isTRUE(actual_rec$use_wild_bootstrap))
  expect_length(actual_rec$alternatives, 0L)
})

test_that("E8-04 Task 3 public-contract oracle freezes local-scope consistency semantics", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_clustering_parity_path(
    "20260324-qa-parity-e8-04-clustering-task23-contract.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing clustering public-contract oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_clustering_oracle(
    "20260324-qa-parity-e8-04-clustering-task23-contract.json"
  )

  hierarchical_case <- oracle$cases$hierarchical_consistency
  hierarchical_fixture <- read_clustering_fixture(hierarchical_case$fixture_csv)
  hierarchical_counts <- compute_cluster_unique_counts(
    hierarchical_fixture,
    cluster_var = hierarchical_case$cluster_var,
    treat_var = hierarchical_case$gvar
  )
  expected_hierarchical_counts <- as_named_int_vector(
    hierarchical_case$cluster_unique_counts
  )

  expect_identical(
    .detect_treatment_variation_level(
      hierarchical_fixture,
      ivar = hierarchical_case$ivar,
      potential_cluster_vars = hierarchical_case$cluster_var,
      gvar = hierarchical_case$gvar,
      d = NULL
    ),
    hierarchical_case$python_consistency$`treatment_variation_level`
  )
  expect_identical(
    .determine_cluster_level(
      hierarchical_fixture,
      hierarchical_case$ivar,
      hierarchical_case$cluster_var
    ),
    hierarchical_case$python_consistency$cluster_level
  )
  expect_equal(
    hierarchical_counts[order(names(hierarchical_counts))],
    expected_hierarchical_counts[order(names(expected_hierarchical_counts))]
  )
  expect_equal(
    sum(hierarchical_counts > 1L) / length(hierarchical_counts) * 100,
    hierarchical_case$python_consistency$pct_clusters_with_variation,
    tolerance = 1e-12
  )

  actual_hierarchical <- check_clustering_consistency(
    hierarchical_fixture,
    ivar = hierarchical_case$ivar,
    cluster_var = hierarchical_case$cluster_var,
    gvar = hierarchical_case$gvar,
    d = NULL,
    verbose = FALSE
  )
  expect_identical(
    actual_hierarchical$is_consistent,
    hierarchical_case$python_consistency$is_consistent
  )
  expect_identical(
    actual_hierarchical$treatment_variation_level,
    hierarchical_case$python_consistency$`treatment_variation_level`
  )
  expect_identical(
    actual_hierarchical$cluster_level,
    hierarchical_case$python_consistency$cluster_level
  )
  expect_identical(
    actual_hierarchical$n_clusters,
    hierarchical_case$python_consistency$n_clusters
  )
  expect_identical(
    actual_hierarchical$n_inconsistent,
    hierarchical_case$python_consistency$n_treatment_changes_within_cluster
  )
  expect_equal(
    actual_hierarchical$pct_inconsistent,
    hierarchical_case$python_consistency$pct_clusters_with_variation,
    tolerance = 1e-12
  )
  expect_identical(
    actual_hierarchical$recommendation,
    hierarchical_case$python_consistency$recommendation
  )
  expect_identical(
    actual_hierarchical$details,
    hierarchical_case$python_consistency$details
  )
  expect_identical(actual_hierarchical$inconsistent_clusters, character(0))

  never_treated_case <- oracle$cases$never_treated_consistency
  never_treated_fixture <- read_clustering_fixture(never_treated_case$fixture_csv)
  never_treated_counts <- compute_cluster_unique_counts(
    never_treated_fixture,
    cluster_var = never_treated_case$cluster_var,
    treat_var = never_treated_case$gvar
  )
  expected_never_treated_counts <- as_named_int_vector(
    never_treated_case$cluster_unique_counts
  )

  expect_identical(
    .detect_treatment_variation_level(
      never_treated_fixture,
      ivar = never_treated_case$ivar,
      potential_cluster_vars = never_treated_case$cluster_var,
      gvar = never_treated_case$gvar,
      d = NULL
    ),
    never_treated_case$python_consistency$`treatment_variation_level`
  )
  expect_identical(
    .determine_cluster_level(
      never_treated_fixture,
      never_treated_case$ivar,
      never_treated_case$cluster_var
    ),
    never_treated_case$python_consistency$cluster_level
  )
  expect_equal(
    never_treated_counts[order(names(never_treated_counts))],
    expected_never_treated_counts[order(names(expected_never_treated_counts))]
  )
  expect_equal(
    sum(never_treated_counts > 1L) / length(never_treated_counts) * 100,
    never_treated_case$python_consistency$pct_clusters_with_variation,
    tolerance = 1e-12
  )
  expect_false(isTRUE(never_treated_case$python_consistency$is_consistent))

  actual_never_treated <- check_clustering_consistency(
    never_treated_fixture,
    ivar = never_treated_case$ivar,
    cluster_var = never_treated_case$cluster_var,
    gvar = never_treated_case$gvar,
    d = NULL,
    verbose = FALSE
  )
  expect_identical(
    actual_never_treated$is_consistent,
    never_treated_case$python_consistency$is_consistent
  )
  expect_identical(
    actual_never_treated$treatment_variation_level,
    never_treated_case$python_consistency$`treatment_variation_level`
  )
  expect_identical(
    actual_never_treated$cluster_level,
    never_treated_case$python_consistency$cluster_level
  )
  expect_identical(
    actual_never_treated$n_clusters,
    never_treated_case$python_consistency$n_clusters
  )
  expect_identical(
    actual_never_treated$n_inconsistent,
    never_treated_case$python_consistency$n_treatment_changes_within_cluster
  )
  expect_equal(
    actual_never_treated$pct_inconsistent,
    never_treated_case$python_consistency$pct_clusters_with_variation,
    tolerance = 1e-12
  )
  expect_identical(
    actual_never_treated$recommendation,
    never_treated_case$python_consistency$recommendation
  )
  expect_identical(
    actual_never_treated$details,
    never_treated_case$python_consistency$details
  )
  expect_identical(actual_never_treated$inconsistent_clusters, "C")
})

test_that("E8-04 Layer 3 parity freezes state-policy public workflow", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_clustering_parity_path(
    "20260324-story-worker-e8-04-layer3-state-policy.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing clustering layer-3 oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_clustering_oracle(
    "20260324-story-worker-e8-04-layer3-state-policy.json"
  )
  case <- oracle$cases$state_policy
  fixture <- read_clustering_fixture(case$fixture_csv)

  actual_diag <- diagnose_clustering(
    fixture,
    ivar = case$ivar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = case$gvar,
    d = NULL,
    verbose = FALSE
  )

  expect_identical(
    actual_diag$recommended_cluster_var,
    case$python_diagnose_clustering$recommended_cluster_var
  )
  expect_identical(
    actual_diag$treatment_variation_level,
    case$python_diagnose_clustering$treatment_variation_level
  )
  expect_identical(
    actual_diag$recommendation_reason,
    case$python_diagnose_clustering$recommendation_reason
  )
  expect_identical(
    actual_diag$warnings,
    as_character_vector(case$python_diagnose_clustering$warnings)
  )

  for (cluster_var in names(case$diagnosis_cluster_structure)) {
    expected_stats <- case$diagnosis_cluster_structure[[cluster_var]]
    actual_stats <- actual_diag$cluster_structure[[cluster_var]]

    expect_identical(actual_stats$n_clusters, expected_stats$n_clusters)
    expect_identical(
      actual_stats$level_relative_to_unit,
      expected_stats$level_relative_to_unit
    )
    expect_identical(
      actual_stats$treatment_varies_within,
      expected_stats$treatment_varies_within
    )
    expect_equal(
      actual_stats$reliability_score,
      expected_stats$reliability_score,
      tolerance = 1e-12
    )
  }

  actual_rec <- recommend_clustering(
    fixture,
    ivar = case$ivar,
    tvar = case$tvar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = case$gvar,
    d = NULL,
    min_clusters = case$min_clusters,
    verbose = FALSE
  )

  expect_identical(
    actual_rec$recommended_var,
    case$python_recommend_clustering$recommended_var
  )
  expect_identical(actual_rec$n_clusters, case$python_recommend_clustering$n_clusters)
  expect_identical(
    actual_rec$n_treated_clusters,
    case$python_recommend_clustering$n_treated_clusters
  )
  expect_identical(
    actual_rec$n_control_clusters,
    case$python_recommend_clustering$n_control_clusters
  )
  expect_equal(
    actual_rec$confidence,
    case$python_recommend_clustering$confidence,
    tolerance = 1e-12
  )
  expect_identical(
    actual_rec$reasons,
    unlist(case$python_recommend_clustering$reasons, use.names = FALSE)
  )
  expect_identical(
    actual_rec$warnings,
    as_character_vector(case$python_recommend_clustering$warnings)
  )
  expect_identical(
    actual_rec$use_wild_bootstrap,
    case$python_recommend_clustering$use_wild_bootstrap
  )
  expect_identical(
    actual_rec$wild_bootstrap_reason,
    case$python_recommend_clustering$wild_bootstrap_reason
  )
  expect_length(
    actual_rec$alternatives,
    length(case$python_recommend_clustering$alternatives)
  )
  for (idx in seq_along(actual_rec$alternatives)) {
    expect_identical(
      actual_rec$alternatives[[idx]]$var,
      case$python_recommend_clustering$alternatives[[idx]]$var
    )
    expect_identical(
      actual_rec$alternatives[[idx]]$n_clusters,
      case$python_recommend_clustering$alternatives[[idx]]$n_clusters
    )
    expect_equal(
      actual_rec$alternatives[[idx]]$reliability_score,
      case$python_recommend_clustering$alternatives[[idx]]$reliability_score,
      tolerance = 1e-12
    )
    expect_identical(
      actual_rec$alternatives[[idx]]$reason,
      case$python_recommend_clustering$alternatives[[idx]]$reason
    )
  }

  actual_consistency <- check_clustering_consistency(
    fixture,
    ivar = case$ivar,
    cluster_var = case$consistency_cluster_var,
    gvar = case$gvar,
    d = NULL,
    verbose = FALSE
  )

  expect_identical(
    actual_consistency$is_consistent,
    case$python_consistency$is_consistent
  )
  expect_identical(
    actual_consistency$treatment_variation_level,
    case$python_consistency$treatment_variation_level
  )
  expect_identical(
    actual_consistency$cluster_level,
    case$python_consistency$cluster_level
  )
  expect_identical(
    actual_consistency$n_clusters,
    case$python_consistency$n_clusters
  )
  expect_identical(
    actual_consistency$n_inconsistent,
    case$python_consistency$n_treatment_changes_within_cluster
  )
  expect_equal(
    actual_consistency$pct_inconsistent,
    case$python_consistency$pct_clusters_with_variation,
    tolerance = 1e-12
  )
  expect_identical(
    actual_consistency$recommendation,
    case$python_consistency$recommendation
  )
  expect_identical(
    actual_consistency$details,
    case$python_consistency$details
  )
})

test_that("E8-04 Layer 3 parity freezes few-cluster public workflow", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_clustering_parity_path(
    "20260324-qa-parity-e8-04-layer3-few-cluster.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing clustering few-cluster oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_clustering_oracle(
    "20260324-qa-parity-e8-04-layer3-few-cluster.json"
  )
  case <- oracle$cases$few_cluster
  fixture <- read_clustering_fixture(case$fixture_csv)

  actual_diag <- diagnose_clustering(
    fixture,
    ivar = case$ivar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = case$gvar,
    d = NULL,
    verbose = FALSE
  )

  expect_identical(
    actual_diag$recommended_cluster_var,
    case$python_diagnose_clustering$recommended_cluster_var
  )
  expect_identical(
    actual_diag$treatment_variation_level,
    case$python_diagnose_clustering$treatment_variation_level
  )
  expect_identical(
    actual_diag$recommendation_reason,
    case$python_diagnose_clustering$recommendation_reason
  )
  expect_identical(
    actual_diag$warnings,
    as_character_vector(case$python_diagnose_clustering$warnings)
  )

  for (cluster_var in names(case$diagnosis_cluster_structure)) {
    expected_stats <- case$diagnosis_cluster_structure[[cluster_var]]
    actual_stats <- actual_diag$cluster_structure[[cluster_var]]

    expect_identical(actual_stats$n_clusters, expected_stats$n_clusters)
    expect_identical(
      actual_stats$level_relative_to_unit,
      expected_stats$level_relative_to_unit
    )
    expect_identical(
      actual_stats$treatment_varies_within,
      expected_stats$treatment_varies_within
    )
    expect_equal(
      actual_stats$reliability_score,
      expected_stats$reliability_score,
      tolerance = 1e-12
    )
  }

  actual_rec <- recommend_clustering(
    fixture,
    ivar = case$ivar,
    tvar = case$tvar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = case$gvar,
    d = NULL,
    min_clusters = case$min_clusters,
    verbose = FALSE
  )

  expect_identical(
    actual_rec$recommended_var,
    case$python_recommend_clustering$recommended_var
  )
  expect_identical(actual_rec$n_clusters, case$python_recommend_clustering$n_clusters)
  expect_identical(
    actual_rec$n_treated_clusters,
    case$python_recommend_clustering$n_treated_clusters
  )
  expect_identical(
    actual_rec$n_control_clusters,
    case$python_recommend_clustering$n_control_clusters
  )
  expect_equal(
    actual_rec$confidence,
    case$python_recommend_clustering$confidence,
    tolerance = 1e-12
  )
  expect_identical(
    actual_rec$reasons,
    unlist(case$python_recommend_clustering$reasons, use.names = FALSE)
  )
  expect_identical(
    actual_rec$warnings,
    as_character_vector(case$python_recommend_clustering$warnings)
  )
  expect_identical(
    actual_rec$use_wild_bootstrap,
    case$python_recommend_clustering$use_wild_bootstrap
  )
  expect_identical(
    actual_rec$wild_bootstrap_reason,
    case$python_recommend_clustering$wild_bootstrap_reason
  )
  expect_length(
    actual_rec$alternatives,
    length(case$python_recommend_clustering$alternatives)
  )
  for (idx in seq_along(actual_rec$alternatives)) {
    expect_identical(
      actual_rec$alternatives[[idx]]$var,
      case$python_recommend_clustering$alternatives[[idx]]$var
    )
    expect_identical(
      actual_rec$alternatives[[idx]]$n_clusters,
      case$python_recommend_clustering$alternatives[[idx]]$n_clusters
    )
    expect_equal(
      actual_rec$alternatives[[idx]]$reliability_score,
      case$python_recommend_clustering$alternatives[[idx]]$reliability_score,
      tolerance = 1e-12
    )
    expect_identical(
      actual_rec$alternatives[[idx]]$reason,
      case$python_recommend_clustering$alternatives[[idx]]$reason
    )
  }

  actual_consistency <- check_clustering_consistency(
    fixture,
    ivar = case$ivar,
    cluster_var = case$consistency_cluster_var,
    gvar = case$gvar,
    d = NULL,
    verbose = FALSE
  )

  expect_identical(
    actual_consistency$is_consistent,
    case$python_consistency$is_consistent
  )
  expect_identical(
    actual_consistency$treatment_variation_level,
    case$python_consistency$treatment_variation_level
  )
  expect_identical(
    actual_consistency$cluster_level,
    case$python_consistency$cluster_level
  )
  expect_identical(
    actual_consistency$n_clusters,
    case$python_consistency$n_clusters
  )
  expect_identical(
    actual_consistency$n_inconsistent,
    case$python_consistency$n_treatment_changes_within_cluster
  )
  expect_equal(
    actual_consistency$pct_inconsistent,
    case$python_consistency$pct_clusters_with_variation,
    tolerance = 1e-12
  )
  expect_identical(
    actual_consistency$recommendation,
    case$python_consistency$recommendation
  )
  expect_identical(
    actual_consistency$details,
    case$python_consistency$details
  )
})

test_that("E8-04 Layer 3 parity freezes hierarchical public workflow", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_clustering_parity_path(
    "20260324-qa-parity-e8-04-layer3-hierarchical.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing clustering hierarchical oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_clustering_oracle(
    "20260324-qa-parity-e8-04-layer3-hierarchical.json"
  )
  case <- oracle$cases$hierarchical
  fixture <- read_clustering_fixture(case$fixture_csv)

  actual_diag <- diagnose_clustering(
    fixture,
    ivar = case$ivar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = case$gvar,
    d = NULL,
    verbose = FALSE
  )

  expect_identical(
    actual_diag$recommended_cluster_var,
    case$python_diagnose_clustering$recommended_cluster_var
  )
  expect_identical(
    actual_diag$treatment_variation_level,
    case$python_diagnose_clustering$treatment_variation_level
  )
  expect_identical(
    actual_diag$recommendation_reason,
    case$python_diagnose_clustering$recommendation_reason
  )
  expect_identical(
    actual_diag$warnings,
    as_character_vector(case$python_diagnose_clustering$warnings)
  )

  for (cluster_var in names(case$diagnosis_cluster_structure)) {
    expected_stats <- case$diagnosis_cluster_structure[[cluster_var]]
    actual_stats <- actual_diag$cluster_structure[[cluster_var]]

    expect_identical(actual_stats$n_clusters, expected_stats$n_clusters)
    expect_identical(
      actual_stats$level_relative_to_unit,
      expected_stats$level_relative_to_unit
    )
    expect_identical(
      actual_stats$treatment_varies_within,
      expected_stats$treatment_varies_within
    )
    expect_equal(
      actual_stats$reliability_score,
      expected_stats$reliability_score,
      tolerance = 1e-12
    )
  }

  actual_rec <- recommend_clustering(
    fixture,
    ivar = case$ivar,
    tvar = case$tvar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = case$gvar,
    d = NULL,
    min_clusters = case$min_clusters,
    verbose = FALSE
  )

  expect_identical(
    actual_rec$recommended_var,
    case$python_recommend_clustering$recommended_var
  )
  expect_identical(actual_rec$n_clusters, case$python_recommend_clustering$n_clusters)
  expect_identical(
    actual_rec$n_treated_clusters,
    case$python_recommend_clustering$n_treated_clusters
  )
  expect_identical(
    actual_rec$n_control_clusters,
    case$python_recommend_clustering$n_control_clusters
  )
  expect_equal(
    actual_rec$confidence,
    case$python_recommend_clustering$confidence,
    tolerance = 1e-12
  )
  expect_identical(
    actual_rec$reasons,
    unlist(case$python_recommend_clustering$reasons, use.names = FALSE)
  )
  expect_identical(
    actual_rec$warnings,
    as_character_vector(case$python_recommend_clustering$warnings)
  )
  expect_identical(
    actual_rec$use_wild_bootstrap,
    case$python_recommend_clustering$use_wild_bootstrap
  )
  expect_identical(
    actual_rec$wild_bootstrap_reason,
    case$python_recommend_clustering$wild_bootstrap_reason
  )
  expect_length(
    actual_rec$alternatives,
    length(case$python_recommend_clustering$alternatives)
  )
  for (idx in seq_along(actual_rec$alternatives)) {
    expect_identical(
      actual_rec$alternatives[[idx]]$var,
      case$python_recommend_clustering$alternatives[[idx]]$var
    )
    expect_identical(
      actual_rec$alternatives[[idx]]$n_clusters,
      case$python_recommend_clustering$alternatives[[idx]]$n_clusters
    )
    expect_equal(
      actual_rec$alternatives[[idx]]$reliability_score,
      case$python_recommend_clustering$alternatives[[idx]]$reliability_score,
      tolerance = 1e-12
    )
    expect_identical(
      actual_rec$alternatives[[idx]]$reason,
      case$python_recommend_clustering$alternatives[[idx]]$reason
    )
  }

  actual_consistency <- check_clustering_consistency(
    fixture,
    ivar = case$ivar,
    cluster_var = case$consistency_cluster_var,
    gvar = case$gvar,
    d = NULL,
    verbose = FALSE
  )

  expect_identical(
    actual_consistency$is_consistent,
    case$python_consistency$is_consistent
  )
  expect_identical(
    actual_consistency$treatment_variation_level,
    case$python_consistency$treatment_variation_level
  )
  expect_identical(
    actual_consistency$cluster_level,
    case$python_consistency$cluster_level
  )
  expect_identical(
    actual_consistency$n_clusters,
    case$python_consistency$n_clusters
  )
  expect_identical(
    actual_consistency$n_inconsistent,
    case$python_consistency$n_treatment_changes_within_cluster
  )
  expect_equal(
    actual_consistency$pct_inconsistent,
    case$python_consistency$pct_clusters_with_variation,
    tolerance = 1e-12
  )
  expect_identical(
    actual_consistency$recommendation,
    case$python_consistency$recommendation
  )
  expect_identical(
    actual_consistency$details,
    case$python_consistency$details
  )
})

test_that("E8-04 Layer 3 parity freezes smoking common-timing public workflow", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_clustering_parity_path(
    "20260324-story-worker-e8-04-layer3-smoking-common-timing.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing clustering smoking common-timing oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_clustering_oracle(
    "20260324-story-worker-e8-04-layer3-smoking-common-timing.json"
  )
  case <- oracle$cases$smoking_common_timing
  utils::data(list = case$r_dataset, package = "lwdid", envir = environment())

  fixture <- get(case$r_dataset, envir = environment())

  actual_diag <- diagnose_clustering(
    fixture,
    ivar = case$ivar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = NULL,
    d = case$d,
    verbose = FALSE
  )

  expect_identical(
    actual_diag$recommended_cluster_var,
    case$python_diagnose_clustering$recommended_cluster_var
  )
  expect_identical(
    actual_diag$treatment_variation_level,
    case$python_diagnose_clustering$treatment_variation_level
  )
  expect_identical(
    actual_diag$recommendation_reason,
    case$python_diagnose_clustering$recommendation_reason
  )
  expect_identical(
    actual_diag$warnings,
    as_character_vector(case$python_diagnose_clustering$warnings)
  )

  for (cluster_var in names(case$diagnosis_cluster_structure)) {
    expected_stats <- case$diagnosis_cluster_structure[[cluster_var]]
    actual_stats <- actual_diag$cluster_structure[[cluster_var]]

    expect_identical(actual_stats$n_clusters, expected_stats$n_clusters)
    expect_identical(
      actual_stats$level_relative_to_unit,
      expected_stats$level_relative_to_unit
    )
    expect_identical(
      actual_stats$treatment_varies_within,
      expected_stats$treatment_varies_within
    )
    expect_equal(
      actual_stats$reliability_score,
      expected_stats$reliability_score,
      tolerance = 1e-12
    )
  }

  actual_rec <- recommend_clustering(
    fixture,
    ivar = case$ivar,
    tvar = case$tvar,
    potential_cluster_vars = unlist(
      case$potential_cluster_vars,
      use.names = FALSE
    ),
    gvar = NULL,
    d = case$d,
    min_clusters = case$min_clusters,
    verbose = FALSE
  )

  expect_identical(
    actual_rec$recommended_var,
    case$python_recommend_clustering$recommended_var
  )
  expect_identical(actual_rec$n_clusters, case$python_recommend_clustering$n_clusters)
  expect_identical(
    actual_rec$n_treated_clusters,
    case$python_recommend_clustering$n_treated_clusters
  )
  expect_identical(
    actual_rec$n_control_clusters,
    case$python_recommend_clustering$n_control_clusters
  )
  expect_equal(
    actual_rec$confidence,
    case$python_recommend_clustering$confidence,
    tolerance = 1e-12
  )
  expect_identical(
    actual_rec$reasons,
    as_character_vector(case$python_recommend_clustering$reasons)
  )
  expect_identical(
    actual_rec$warnings,
    as_character_vector(case$python_recommend_clustering$warnings)
  )
  expect_identical(
    actual_rec$use_wild_bootstrap,
    case$python_recommend_clustering$use_wild_bootstrap
  )
  expect_identical(
    actual_rec$wild_bootstrap_reason,
    case$python_recommend_clustering$wild_bootstrap_reason
  )
  expect_length(
    actual_rec$alternatives,
    length(case$python_recommend_clustering$alternatives)
  )

  actual_consistency <- check_clustering_consistency(
    fixture,
    ivar = case$ivar,
    cluster_var = case$consistency_cluster_var,
    gvar = NULL,
    d = case$d,
    verbose = FALSE
  )

  expect_identical(
    actual_consistency$is_consistent,
    case$python_consistency$is_consistent
  )
  expect_identical(
    actual_consistency$treatment_variation_level,
    case$python_consistency$treatment_variation_level
  )
  expect_identical(
    actual_consistency$cluster_level,
    case$python_consistency$cluster_level
  )
  expect_identical(
    actual_consistency$n_clusters,
    case$python_consistency$n_clusters
  )
  expect_identical(
    actual_consistency$n_inconsistent,
    case$python_consistency$n_treatment_changes_within_cluster
  )
  expect_equal(
    actual_consistency$pct_inconsistent,
    case$python_consistency$pct_clusters_with_variation,
    tolerance = 1e-12
  )
  expect_identical(
    actual_consistency$recommendation,
    case$python_consistency$recommendation
  )
  expect_identical(
    actual_consistency$details,
    case$python_consistency$details
  )
})

read_clustering_layer4_contract <- function() {
  contract_path <- resolve_clustering_parity_path(
    "e8_04_clustering_layer4_monte_carlo_contract.yaml"
  )

  yaml::yaml.load(
    paste(
      readLines(contract_path, encoding = "UTF-8", warn = FALSE),
      collapse = "\n"
    )
  )
}

extract_clustering_layer4_scenario <- function(contract, scenario_id) {
  scenario_ids <- vapply(
    contract[["scenarios"]],
    function(scenario) scenario[["id"]],
    character(1)
  )
  match_idx <- match(scenario_id, scenario_ids)
  if (is.na(match_idx)) {
    stop(sprintf("Layer 4 scenario '%s' not found in contract", scenario_id))
  }

  contract[["scenarios"]][[match_idx]]
}

expect_metric_in_acceptance_band <- function(metric_value, acceptance) {
  if (!is.null(acceptance$lower)) {
    expect_gte(metric_value, acceptance$lower)
  }
  if (!is.null(acceptance$lower_exclusive)) {
    expect_gt(metric_value, acceptance$lower_exclusive)
  }
  if (!is.null(acceptance$upper)) {
    expect_lte(metric_value, acceptance$upper)
  }
  if (!is.null(acceptance$upper_exclusive)) {
    expect_lt(metric_value, acceptance$upper_exclusive)
  }
}

test_that("E8-04 Layer 4 parity freezes standard coverage_large_g acceptance band", {
  skip_if_not_installed("yaml")

  contract_path <- resolve_clustering_parity_path(
    "e8_04_clustering_layer4_monte_carlo_contract.yaml"
  )
  expect_true(
    file.exists(contract_path),
    info = paste("missing clustering Layer 4 contract:", contract_path)
  )
  if (!file.exists(contract_path)) {
    return(invisible(NULL))
  }

  contract <- read_clustering_layer4_contract()
  scenario <- extract_clustering_layer4_scenario(contract, "coverage_large_g")

  summary <- .run_clustering_monte_carlo_scenario(
    scenario,
    shared_dgp = contract[["shared_dgp"]]
  )

  expect_identical(summary$scenario_id, "coverage_large_g")
  expect_identical(summary$n_simulations, scenario$n_simulations)
  expect_identical(summary$metric_name, "coverage_rate")
  expect_true(is.finite(summary$metric_value))
  expect_metric_in_acceptance_band(
    summary$metric_value,
    scenario$acceptance$coverage_rate
  )
})

test_that("E8-04 Layer 4 parity freezes remaining standard acceptance bands", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("yaml")

  audit_path <- resolve_clustering_parity_path(
    "20260324-qa-parity-e8-04-layer4-standard-scenarios.json"
  )
  expect_true(
    file.exists(audit_path),
    info = paste("missing clustering Layer 4 standard audit:", audit_path)
  )
  if (!file.exists(audit_path)) {
    return(invisible(NULL))
  }

  audit <- read_clustering_oracle(
    "20260324-qa-parity-e8-04-layer4-standard-scenarios.json"
  )
  contract <- read_clustering_layer4_contract()

  expected_ids <- c(
    "coverage_medium_g",
    "coverage_small_g",
    "size_under_null",
    "power_under_alternative",
    "unbalanced_treatment_coverage",
    "heterogeneous_cluster_sizes_coverage"
  )

  audit_ids <- vapply(
    audit$scenarios,
    function(entry) entry$scenario_id,
    character(1)
  )
  expect_identical(audit_ids, expected_ids)

  for (scenario_id in expected_ids) {
    audit_idx <- match(scenario_id, audit_ids)
    audit_entry <- audit$scenarios[[audit_idx]]
    scenario <- extract_clustering_layer4_scenario(contract, scenario_id)

    summary <- .run_clustering_monte_carlo_scenario(
      scenario,
      shared_dgp = contract[["shared_dgp"]]
    )

    expect_identical(summary$scenario_id, audit_entry$scenario_id)
    expect_identical(summary$n_simulations, audit_entry$n_simulations)
    expect_identical(summary$metric_name, audit_entry$metric_name)
    expect_equal(summary$metric_value, audit_entry$metric_value, tolerance = 0)
    expect_metric_in_acceptance_band(
      summary$metric_value,
      scenario$acceptance[[summary$metric_name]]
    )
  }
})

test_that("E8-04 Layer 4 parity freezes wild bootstrap coverage_small_g contract", {
  skip_if_not_installed("yaml")

  contract_path <- resolve_clustering_parity_path(
    "e8_04_clustering_layer4_monte_carlo_contract.yaml"
  )
  expect_true(
    file.exists(contract_path),
    info = paste("missing clustering Layer 4 contract:", contract_path)
  )
  if (!file.exists(contract_path)) {
    return(invisible(NULL))
  }

  contract <- read_clustering_layer4_contract()
  scenario <- extract_clustering_layer4_scenario(
    contract,
    "wild_bootstrap_coverage_small_g"
  )

  summary <- .run_clustering_monte_carlo_scenario(
    scenario,
    shared_dgp = contract[["shared_dgp"]]
  )

  expect_identical(summary$scenario_id, "wild_bootstrap_coverage_small_g")
  expect_identical(summary$n_simulations, scenario$n_simulations)
  expect_identical(summary$metric_name, "coverage_rate")
  expect_identical(summary$requested_n_bootstrap, scenario$requested_n_bootstrap)
  expect_identical(summary$actual_n_bootstrap, scenario$actual_n_bootstrap)
  expect_identical(summary$weight_type, scenario$weight_type)
  expect_identical(summary$impose_null, scenario$impose_null)
  expect_identical(summary$ci_method, scenario$ci_method)
  expect_true(isTRUE(summary$full_enumeration))
  expect_true(is.finite(summary$metric_value))
  expect_metric_in_acceptance_band(
    summary$metric_value,
    scenario$acceptance$coverage_rate
  )
})

test_that("E8-04 Layer 4 parity freezes wild bootstrap relative-gap contract", {
  skip_if_not_installed("yaml")

  contract_path <- resolve_clustering_parity_path(
    "e8_04_clustering_layer4_monte_carlo_contract.yaml"
  )
  expect_true(
    file.exists(contract_path),
    info = paste("missing clustering Layer 4 contract:", contract_path)
  )
  if (!file.exists(contract_path)) {
    return(invisible(NULL))
  }

  contract <- read_clustering_layer4_contract()
  scenario <- extract_clustering_layer4_scenario(
    contract,
    "wild_bootstrap_relative_gap"
  )

  summary <- .run_clustering_monte_carlo_scenario(
    scenario,
    shared_dgp = contract[["shared_dgp"]]
  )

  expect_identical(summary$scenario_id, "wild_bootstrap_relative_gap")
  expect_identical(summary$n_simulations, scenario$n_simulations)
  expect_identical(summary$metric_name, "coverage_rate_gap")
  expect_identical(
    summary$requested_n_bootstrap,
    scenario$estimators$wild$requested_n_bootstrap
  )
  expect_identical(
    summary$actual_n_bootstrap,
    scenario$estimators$wild$actual_n_bootstrap
  )
  expect_identical(summary$weight_type, scenario$estimators$wild$weight_type)
  expect_identical(summary$impose_null, scenario$estimators$wild$impose_null)
  expect_identical(summary$ci_method, scenario$estimators$wild$ci_method)
  expect_true(isTRUE(summary$full_enumeration))
  expect_true(is.list(summary$metrics))
  expect_identical(
    sort(names(summary$metrics)),
    sort(c("coverage_rate_gap", "coverage_std_rate", "coverage_wild_rate"))
  )
  expect_equal(
    summary$metric_value,
    summary$metrics$coverage_wild_rate - summary$metrics$coverage_std_rate,
    tolerance = 0
  )
  expect_metric_in_acceptance_band(
    summary$metric_value,
    scenario$acceptance$min_relative_difference
  )
})
