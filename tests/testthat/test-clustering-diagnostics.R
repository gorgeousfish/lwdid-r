library(testthat)

make_cluster_panel <- function() {
  data.frame(
    id = rep(1:6, each = 5),
    time = rep(1:5, 6),
    gvar = rep(c(3, 3, NA, 0, Inf, 4), each = 5),
    cluster = rep(c("A", "A", "B", "B", "C", "C"), each = 5),
    stringsAsFactors = FALSE
  )
}

make_cluster_hierarchy_panel <- function(
    n_regions = 4L,
    states_per_region = 6L,
    units_per_state = 2L,
    periods = 2L
) {
  state_ids <- seq_len(n_regions * states_per_region)
  unit_ids <- seq_len(length(state_ids) * units_per_state)

  state_lookup <- data.frame(
    id = unit_ids,
    state = rep(sprintf("s%02d", state_ids), each = units_per_state),
    region = rep(
      sprintf("r%02d", seq_len(n_regions)),
      each = states_per_region * units_per_state
    ),
    gvar = rep(
      rep(rep(c(2005, 0), length.out = states_per_region), times = n_regions),
      each = units_per_state
    ),
    stringsAsFactors = FALSE
  )

  panel <- merge(
    expand.grid(
      id = unit_ids,
      time = seq_len(periods),
      KEEP.OUT.ATTRS = FALSE
    ),
    state_lookup,
    by = "id",
    sort = FALSE
  )

  panel[order(panel$id, panel$time), ]
}

make_consistency_panel <- function(n_clusters, inconsistent_clusters = integer(0)) {
  panel <- data.frame(
    id = seq_len(n_clusters * 2L),
    time = 1L,
    gvar = rep(5, n_clusters * 2L),
    cluster = rep(sprintf("c%02d", seq_len(n_clusters)), each = 2L),
    stringsAsFactors = FALSE
  )

  if (length(inconsistent_clusters) > 0L) {
    for (idx in inconsistent_clusters) {
      rows <- ((idx - 1L) * 2L + 1L):((idx - 1L) * 2L + 2L)
      panel$gvar[rows] <- c(5, 6)
    }
  }

  panel
}

test_that("E8-04.1.1: .validate_clustering_inputs accepts valid staggered inputs", {
  panel <- make_cluster_panel()

  expect_invisible(
    .validate_clustering_inputs(
      panel,
      ivar = "id",
      potential_cluster_vars = c("cluster"),
      gvar = "gvar",
      d = NULL
    )
  )
})

test_that("E8-04.1.1: .validate_clustering_inputs rejects invalid inputs", {
  panel <- make_cluster_panel()

  expect_error(
    .validate_clustering_inputs(
      panel,
      ivar = "id",
      potential_cluster_vars = character(0),
      gvar = "gvar",
      d = NULL
    ),
    "At least one potential cluster variable must be specified"
  )

  expect_error(
    .validate_clustering_inputs(
      panel,
      ivar = "id",
      potential_cluster_vars = c("cluster"),
      gvar = NULL,
      d = NULL
    ),
    "Either gvar or d must be specified"
  )

  expect_error(
    .validate_clustering_inputs(
      panel,
      ivar = "missing_id",
      potential_cluster_vars = c("cluster"),
      gvar = "gvar",
      d = NULL
    ),
    "Unit variable 'missing_id' not found in data"
  )
})

test_that("E8-04.1.3: .determine_cluster_level distinguishes higher same and lower", {
  same_level <- data.frame(
    id = rep(1:3, each = 2),
    cluster_same = rep(1:3, each = 2),
    stringsAsFactors = FALSE
  )

  higher_level <- data.frame(
    id = rep(1:6, each = 2),
    state = rep(c("A", "A", "B", "B", "C", "C"), each = 2),
    stringsAsFactors = FALSE
  )

  lower_level <- data.frame(
    id = c(1, 1, 2, 2, 3, 3),
    subcluster = c("a1", "a2", "b1", "b2", "c1", "c2"),
    stringsAsFactors = FALSE
  )

  expect_equal(
    .determine_cluster_level(higher_level, "id", "state"),
    "higher"
  )
  expect_equal(
    .determine_cluster_level(same_level, "id", "cluster_same"),
    "same"
  )
  expect_equal(
    .determine_cluster_level(lower_level, "id", "subcluster"),
    "lower"
  )
})

test_that("E8-04.1.2: .analyze_cluster_var uses observation counts and never-treated rules", {
  panel <- make_cluster_panel()

  result <- .analyze_cluster_var(
    panel,
    ivar = "id",
    cluster_var = "cluster",
    gvar = "gvar",
    d = NULL
  )

  expect_equal(result$var_name, "cluster")
  expect_equal(result$n_clusters, 3L)
  expect_equal(result$cluster_sizes, c(A = 10L, B = 10L, C = 10L))
  expect_equal(result$min_size, 10L)
  expect_equal(result$max_size, 10L)
  expect_equal(result$mean_size, 10)
  expect_equal(result$median_size, 10)
  expect_equal(result$cv, 0, tolerance = 1e-12)
  expect_equal(result$n_treated_clusters, 2L)
  expect_equal(result$n_control_clusters, 1L)
  expect_true(result$treatment_varies_within)
  expect_equal(result$n_clusters_with_variation, 1L)
  expect_false(result$is_nested_in_unit)
  expect_equal(result$level_relative_to_unit, "higher")
  expect_equal(result$units_per_cluster, 2)
  expect_equal(result$balance_score, 2 / 3, tolerance = 1e-12)
  expect_equal(result$cv_score, 1, tolerance = 1e-12)
  expect_equal(result$reliability_score, 0.43, tolerance = 1e-12)
  expect_true(result$is_valid)
})

test_that("E8-04.1.2: .analyze_cluster_var counts staggered cohort variation within clusters", {
  staggered_panel <- data.frame(
    id = rep(1:4, each = 4),
    time = rep(1:4, 4),
    gvar = rep(c(2005, 2010, 2005, 2005), each = 4),
    cluster = rep(c("A", "A", "B", "B"), each = 4),
    stringsAsFactors = FALSE
  )

  result <- .analyze_cluster_var(
    staggered_panel,
    ivar = "id",
    cluster_var = "cluster",
    gvar = "gvar",
    d = NULL
  )

  expect_true(result$treatment_varies_within)
  expect_equal(result$n_clusters_with_variation, 1L)
})

test_that("E8-04.2.1: .detect_treatment_variation_level finds first non-varying level and falls back to ivar", {
  staggered_panel <- make_cluster_hierarchy_panel()

  expect_equal(
    .detect_treatment_variation_level(
      staggered_panel,
      ivar = "id",
      potential_cluster_vars = c("state", "region"),
      gvar = "gvar",
      d = NULL
    ),
    "state"
  )

  common_timing_panel <- transform(
    staggered_panel,
    d = as.integer(time >= 2L & gvar == 2005)
  )

  expect_equal(
    .detect_treatment_variation_level(
      common_timing_panel,
      ivar = "id",
      potential_cluster_vars = c("region", "state"),
      gvar = NULL,
      d = "d"
    ),
    "id"
  )
})

test_that("E8-04.2.2: .generate_clustering_recommendation uses treatment-level immediate return and dual-key fallback", {
  valid_state <- list(
    var_name = "state",
    n_clusters = 24L,
    reliability_score = 0.62,
    is_valid = TRUE,
    treatment_varies_within = FALSE,
    n_clusters_with_variation = 0L
  )
  strong_county <- list(
    var_name = "county",
    n_clusters = 48L,
    reliability_score = 0.95,
    is_valid = TRUE,
    treatment_varies_within = FALSE,
    n_clusters_with_variation = 0L
  )

  immediate <- .generate_clustering_recommendation(
    cluster_structure = list(state = valid_state, county = strong_county),
    treatment_level = "state"
  )
  expect_identical(immediate$recommended_var, "state")
  expect_identical(immediate$warnings, character(0))

  low_state <- list(
    var_name = "state",
    n_clusters = 15L,
    reliability_score = 0.92,
    is_valid = TRUE,
    treatment_varies_within = FALSE,
    n_clusters_with_variation = 0L
  )
  adequate_county <- list(
    var_name = "county",
    n_clusters = 25L,
    reliability_score = 0.55,
    is_valid = TRUE,
    treatment_varies_within = FALSE,
    n_clusters_with_variation = 0L
  )

  fallback <- .generate_clustering_recommendation(
    cluster_structure = list(state = low_state, county = adequate_county),
    treatment_level = "state"
  )
  expect_identical(fallback$recommended_var, "county")
  expect_true(any(grepl("state", fallback$warnings, fixed = TRUE)))

  all_small <- .generate_clustering_recommendation(
    cluster_structure = list(state = low_state),
    treatment_level = "state"
  )
  expect_identical(all_small$recommended_var, "state")
  expect_identical(
    all_small$warnings,
    c(
      "Treatment varies at state level but only 15 clusters available.",
      paste0(
        "Recommended clustering has only 15 clusters. ",
        "Consider wild cluster bootstrap for reliable inference."
      )
    )
  )
})

test_that("E8-04.3.1: diagnose_clustering returns the public diagnosis contract", {
  panel <- make_cluster_hierarchy_panel()

  expect_silent(
    diag <- diagnose_clustering(
      panel,
      ivar = "id",
      potential_cluster_vars = c("region", "state"),
      gvar = "gvar",
      verbose = FALSE
    )
  )

  expect_s3_class(diag, "lwdid_clustering_diagnosis")
  expect_named(
    diag,
    c(
      "cluster_structure",
      "recommended_cluster_var",
      "recommendation_reason",
      "treatment_variation_level",
      "warnings"
    )
  )
  expect_identical(diag$treatment_variation_level, "state")
  expect_identical(diag$recommended_cluster_var, "state")
  expect_type(diag$cluster_structure, "list")
  expect_true(all(c("region", "state") %in% names(diag$cluster_structure)))

  expect_output(
    diagnose_clustering(
      panel,
      ivar = "id",
      potential_cluster_vars = c("region", "state"),
      gvar = "gvar",
      verbose = TRUE
    ),
    "Treatment variation level"
  )
})

test_that("E8-04.3.2: recommend_clustering returns detailed public recommendation contract", {
  panel <- make_cluster_hierarchy_panel(
    n_regions = 1L,
    states_per_region = 5L,
    units_per_state = 4L,
    periods = 3L
  )

  expect_silent(
    rec <- recommend_clustering(
      panel,
      ivar = "id",
      tvar = "time",
      potential_cluster_vars = c("state"),
      gvar = "gvar",
      verbose = FALSE
    )
  )

  expect_s3_class(rec, "lwdid_clustering_recommendation")
  expect_named(
    rec,
    c(
      "recommended_var",
      "n_clusters",
      "n_treated_clusters",
      "n_control_clusters",
      "confidence",
      "reasons",
      "alternatives",
      "warnings",
      "use_wild_bootstrap",
      "wild_bootstrap_reason"
    )
  )
  expect_identical(rec$recommended_var, "state")
  expect_identical(rec$n_clusters, 5L)
  expect_identical(rec$n_treated_clusters, 3L)
  expect_identical(rec$n_control_clusters, 2L)
  expect_equal(rec$confidence, 0.49, tolerance = 1e-12)
  expect_true(isTRUE(rec$use_wild_bootstrap))
  expect_match(rec$wild_bootstrap_reason, "Wild cluster bootstrap", fixed = TRUE)
  expect_length(rec$alternatives, 0L)
  expect_true(any(grepl("Treatment varies at state level", rec$reasons, fixed = TRUE)))
  expect_true(any(grepl("Limited clusters (5)", rec$reasons, fixed = TRUE)))
  expect_true(any(grepl("Good balance between treated (3) and control (2) clusters", rec$reasons, fixed = TRUE)))
  expect_true(any(grepl("Only 5 clusters", rec$warnings, fixed = TRUE)))

  expect_output(print(rec), "Clustering Recommendation")
  expect_true("recommend_clustering" %in% getNamespaceExports("lwdid"))
})

test_that("E8-04.3.3: check_clustering_consistency keeps local scope and dropna semantics", {
  panel <- make_cluster_panel()

  expect_silent(
    consistency <- check_clustering_consistency(
      panel,
      ivar = "id",
      cluster_var = "cluster",
      gvar = "gvar",
      verbose = FALSE
    )
  )

  expect_named(
    consistency,
    c(
      "is_consistent",
      "treatment_variation_level",
      "cluster_level",
      "n_clusters",
      "n_inconsistent",
      "pct_inconsistent",
      "inconsistent_clusters",
      "recommendation",
      "details"
    )
  )
  expect_false(isTRUE(consistency$is_consistent))
  expect_identical(consistency$treatment_variation_level, "id")
  expect_identical(consistency$cluster_level, "higher")
  expect_identical(consistency$n_clusters, 3L)
  expect_identical(consistency$n_inconsistent, 1L)
  expect_equal(consistency$pct_inconsistent, 100 / 3, tolerance = 1e-12)
  expect_identical(consistency$inconsistent_clusters, "C")
  expect_match(consistency$recommendation, "33.3% of clusters", fixed = TRUE)
  expect_match(consistency$details, "Treatment variation level: id", fixed = TRUE)

  expect_output(print(diagnose_clustering(
    panel,
    ivar = "id",
    potential_cluster_vars = c("cluster"),
    gvar = "gvar",
    verbose = FALSE
  )), "Clustering Diagnostics")
  expect_true("check_clustering_consistency" %in% getNamespaceExports("lwdid"))
})

test_that("E8-04.4.4: consistency threshold uses strict less-than-five percent boundary", {
  four_percent <- check_clustering_consistency(
    make_consistency_panel(25L, 1L),
    ivar = "id",
    cluster_var = "cluster",
    gvar = "gvar",
    verbose = FALSE
  )
  expect_equal(four_percent$pct_inconsistent, 4, tolerance = 1e-12)
  expect_true(isTRUE(four_percent$is_consistent))

  five_percent <- check_clustering_consistency(
    make_consistency_panel(20L, 1L),
    ivar = "id",
    cluster_var = "cluster",
    gvar = "gvar",
    verbose = FALSE
  )
  expect_equal(five_percent$pct_inconsistent, 5, tolerance = 1e-12)
  expect_false(isTRUE(five_percent$is_consistent))

  zero_percent <- check_clustering_consistency(
    make_consistency_panel(20L),
    ivar = "id",
    cluster_var = "cluster",
    gvar = "gvar",
    verbose = FALSE
  )
  expect_equal(zero_percent$pct_inconsistent, 0, tolerance = 1e-12)
  expect_true(isTRUE(zero_percent$is_consistent))

  all_inconsistent <- check_clustering_consistency(
    make_consistency_panel(2L, c(1L, 2L)),
    ivar = "id",
    cluster_var = "cluster",
    gvar = "gvar",
    verbose = FALSE
  )
  expect_equal(all_inconsistent$pct_inconsistent, 100, tolerance = 1e-12)
  expect_false(isTRUE(all_inconsistent$is_consistent))
})

test_that("E8-04.4.12 and 5.2: null recommendations print '无' for invalid clustering choices", {
  invalid_panel <- data.frame(
    id = c(1, 1, 2, 2, 3, 3),
    time = c(1, 2, 1, 2, 1, 2),
    gvar = c(3, 3, Inf, Inf, 4, 4),
    nested_var = c("a1", "a2", "b1", "b2", "c1", "c2"),
    stringsAsFactors = FALSE
  )

  diag <- diagnose_clustering(
    invalid_panel,
    ivar = "id",
    potential_cluster_vars = c("nested_var"),
    gvar = "gvar",
    verbose = FALSE
  )
  expect_null(diag$recommended_cluster_var)
  expect_output(
    print(diag),
    "Recommended cluster variable: None"
  )

  rec <- recommend_clustering(
    invalid_panel,
    ivar = "id",
    tvar = "time",
    potential_cluster_vars = c("nested_var"),
    gvar = "gvar",
    verbose = FALSE
  )
  expect_null(rec$recommended_var)
  expect_identical(rec$confidence, 0)
  expect_output(
    print(rec),
    "Recommended cluster variable: None"
  )
})

test_that("E8-04.5.2: boundary edge cases stay executable", {
  single_candidate <- data.frame(
    id = rep(1:12, each = 2),
    time = rep(1:2, 12),
    gvar = rep(c(rep(5, 6), rep(Inf, 6)), each = 2),
    region = rep(c(rep("A", 6), rep("B", 6)), each = 2),
    stringsAsFactors = FALSE
  )

  result_single <- diagnose_clustering(
    single_candidate,
    ivar = "id",
    potential_cluster_vars = c("region"),
    gvar = "gvar",
    verbose = FALSE
  )
  expect_s3_class(result_single, "lwdid_clustering_diagnosis")
  expect_length(result_single$cluster_structure, 1L)
  expect_identical(names(result_single$cluster_structure), "region")
  expect_identical(result_single$recommended_cluster_var, "region")

  expect_error(
    diagnose_clustering(
      single_candidate,
      ivar = "id",
      potential_cluster_vars = character(0),
      gvar = "gvar",
      verbose = FALSE
    ),
    "At least one potential cluster variable must be specified"
  )
  expect_error(
    diagnose_clustering(
      single_candidate,
      ivar = "id",
      potential_cluster_vars = "missing_cluster",
      gvar = "gvar",
      verbose = FALSE
    ),
    "Cluster variable 'missing_cluster' not found in data"
  )

  one_cluster_data <- single_candidate
  one_cluster_data$single_cluster <- "A"
  result_g1 <- .analyze_cluster_var(
    one_cluster_data,
    ivar = "id",
    cluster_var = "single_cluster",
    gvar = "gvar",
    d = NULL
  )
  expect_identical(result_g1$n_clusters, 1L)
  expect_false(isTRUE(result_g1$is_valid))

  same_level_data <- single_candidate
  same_level_data$unit_cluster <- same_level_data$id
  result_same <- .analyze_cluster_var(
    same_level_data,
    ivar = "id",
    cluster_var = "unit_cluster",
    gvar = "gvar",
    d = NULL
  )
  expect_identical(result_same$level_relative_to_unit, "same")

  equal_data <- data.frame(
    id = rep(1:10, each = 5),
    time = rep(1:5, 10),
    gvar = rep(c(rep(5, 5), rep(Inf, 5)), each = 5),
    equal_cluster = rep(c("A", "A", "B", "B", "C", "C", "D", "D", "E", "E"), each = 5),
    stringsAsFactors = FALSE
  )
  result_equal <- .analyze_cluster_var(
    equal_data,
    ivar = "id",
    cluster_var = "equal_cluster",
    gvar = "gvar",
    d = NULL
  )
  expect_equal(result_equal$cv, 0, tolerance = 1e-12)
  expect_equal(result_equal$cv_score, 1, tolerance = 1e-12)

  all_nt_data <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(1:3, 4),
    gvar = rep(Inf, 12),
    cluster = rep(c("A", "A", "B", "B"), each = 3),
    stringsAsFactors = FALSE
  )
  result_all_nt <- .analyze_cluster_var(
    all_nt_data,
    ivar = "id",
    cluster_var = "cluster",
    gvar = "gvar",
    d = NULL
  )
  expect_identical(result_all_nt$n_treated_clusters, 0L)
  expect_equal(result_all_nt$balance_score, 0, tolerance = 1e-12)
})

test_that("E8-04.4.1, 5.1, and 5.4: reliability score math stays exact", {
  calc_score <- function(G, G_treated, G_control, cv) {
    g_score <- min(G / 50, 1.0)
    balance_score <- if (G > 0L) {
      min(min(G_treated, G_control) / (G / 2), 1.0)
    } else {
      0
    }
    cv_score <- max(0, 1 - cv / 2)

    0.5 * g_score + 0.3 * balance_score + 0.2 * cv_score
  }

  expect_equal(calc_score(40, 20, 20, 0.5), 0.85, tolerance = 1e-10)
  expect_equal(calc_score(10, 8, 2, 1.5), 0.27, tolerance = 1e-10)
  expect_equal(calc_score(100, 50, 50, 0.1), 0.99, tolerance = 1e-10)
  expect_equal(calc_score(1, 0, 1, 0), 0.21, tolerance = 1e-10)
  expect_equal(calc_score(50, 50, 0, 0), 0.70, tolerance = 1e-10)
  expect_equal(calc_score(20, 10, 10, 3.0), 0.50, tolerance = 1e-10)

  balance_cases <- list(
    c(G = 40, treated = 20),
    c(G = 10, treated = 8),
    c(G = 100, treated = 50),
    c(G = 20, treated = 1),
    c(G = 6, treated = 3),
    c(G = 5, treated = 0)
  )

  for (case in balance_cases) {
    G <- unname(case[["G"]])
    treated <- unname(case[["treated"]])
    control <- G - treated
    formula_one <- min(min(treated, control) / (G / 2), 1.0)
    formula_two <- max(0, 1 - abs(treated / G - 0.5) * 2)

    expect_equal(
      formula_one,
      formula_two,
      tolerance = 1e-10,
      info = sprintf("G=%d treated=%d", G, treated)
    )
  }

  panel <- make_cluster_panel()
  result <- .analyze_cluster_var(
    panel,
    ivar = "id",
    cluster_var = "cluster",
    gvar = "gvar",
    d = NULL
  )

  expect_equal(
    result$reliability_score,
    calc_score(
      result$n_clusters,
      result$n_treated_clusters,
      result$n_control_clusters,
      result$cv
    ),
    tolerance = 1e-12
  )
})

test_that("E8-04.4.8: is_valid respects nestedness cluster count and level", {
  nested_panel <- data.frame(
    id = c(1, 1, 2, 2, 3, 3),
    time = c(1, 2, 1, 2, 1, 2),
    gvar = c(5, 5, Inf, Inf, 5, 5),
    nested_var = c("a1", "a2", "b1", "b2", "c1", "c2"),
    stringsAsFactors = FALSE
  )
  nested_result <- .analyze_cluster_var(
    nested_panel,
    ivar = "id",
    cluster_var = "nested_var",
    gvar = "gvar",
    d = NULL
  )
  expect_true(isTRUE(nested_result$is_nested_in_unit))
  expect_identical(nested_result$level_relative_to_unit, "lower")
  expect_false(isTRUE(nested_result$is_valid))

  one_cluster_panel <- data.frame(
    id = rep(1:4, each = 2),
    time = rep(1:2, 4),
    gvar = rep(c(5, 5, Inf, Inf), each = 2),
    single_cluster = "A",
    stringsAsFactors = FALSE
  )
  one_cluster_result <- .analyze_cluster_var(
    one_cluster_panel,
    ivar = "id",
    cluster_var = "single_cluster",
    gvar = "gvar",
    d = NULL
  )
  expect_identical(one_cluster_result$n_clusters, 1L)
  expect_false(isTRUE(one_cluster_result$is_valid))

  valid_panel <- data.frame(
    id = rep(1:4, each = 2),
    time = rep(1:2, 4),
    gvar = rep(c(5, 5, Inf, Inf), each = 2),
    unit_cluster = rep(1:4, each = 2),
    stringsAsFactors = FALSE
  )
  valid_result <- .analyze_cluster_var(
    valid_panel,
    ivar = "id",
    cluster_var = "unit_cluster",
    gvar = "gvar",
    d = NULL
  )
  expect_false(isTRUE(valid_result$is_nested_in_unit))
  expect_identical(valid_result$level_relative_to_unit, "same")
  expect_true(isTRUE(valid_result$is_valid))
})

test_that("E8-04.4.3, 4.7, and 4.9: recommendation boundaries stay aligned", {
  dual_key_choice <- .generate_clustering_recommendation(
    cluster_structure = list(
      state20 = list(
        var_name = "state20",
        n_clusters = 20L,
        reliability_score = 0.41,
        is_valid = TRUE,
        treatment_varies_within = FALSE,
        n_clusters_with_variation = 0L
      ),
      county19 = list(
        var_name = "county19",
        n_clusters = 19L,
        reliability_score = 0.95,
        is_valid = TRUE,
        treatment_varies_within = FALSE,
        n_clusters_with_variation = 0L
      )
    ),
    treatment_level = "region"
  )
  expect_identical(dual_key_choice$recommended_var, "state20")
  expect_false(any(grepl(
    "wild cluster bootstrap",
    dual_key_choice$warnings,
    ignore.case = TRUE
  )))

  small_cluster_choice <- .generate_clustering_recommendation(
    cluster_structure = list(
      state = list(
        var_name = "state",
        n_clusters = 19L,
        reliability_score = 0.72,
        is_valid = TRUE,
        treatment_varies_within = FALSE,
        n_clusters_with_variation = 0L
      )
    ),
    treatment_level = "state"
  )
  expect_identical(small_cluster_choice$recommended_var, "state")
  expect_identical(
    small_cluster_choice$warnings,
    c(
      "Treatment varies at state level but only 19 clusters available.",
      paste0(
        "Recommended clustering has only 19 clusters. ",
        "Consider wild cluster bootstrap for reliable inference."
      )
    )
  )

  adequate_cluster_choice <- .generate_clustering_recommendation(
    cluster_structure = list(
      state = list(
        var_name = "state",
        n_clusters = 20L,
        reliability_score = 0.42,
        is_valid = TRUE,
        treatment_varies_within = FALSE,
        n_clusters_with_variation = 0L
      )
    ),
    treatment_level = "state"
  )
  expect_identical(adequate_cluster_choice$recommended_var, "state")
  expect_false(any(grepl(
    "wild cluster bootstrap",
    adequate_cluster_choice$warnings,
    ignore.case = TRUE
  )))

  panel <- make_cluster_hierarchy_panel(
    n_regions = 2L,
    states_per_region = 12L,
    units_per_state = 2L,
    periods = 2L
  )
  state_index <- as.integer(sub("^s", "", panel$state))
  panel$division <- sprintf("d%02d", ceiling(state_index / 6))
  panel$zone <- sprintf("z%02d", ceiling(state_index / 4))
  panel$block <- sprintf("b%02d", ceiling(state_index / 3))

  rec <- recommend_clustering(
    panel,
    ivar = "id",
    tvar = "time",
    potential_cluster_vars = c("state", "block", "zone", "division"),
    gvar = "gvar",
    verbose = FALSE
  )

  expect_identical(rec$recommended_var, "state")
  expect_length(rec$alternatives, 2L)
})
