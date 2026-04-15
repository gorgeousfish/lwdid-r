library(testthat)

resolve_e8_08_package_dir <- function() {
  candidates <- c(
    ".",
    "lwdid-r",
    file.path("..", ".."),
    file.path("..", "..", "..", "lwdid-r")
  )
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  existing <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]

  if (length(existing) > 0L) {
    return(existing[[1L]])
  }

  NULL
}

if (!exists("lwdid_sensitivity", mode = "function")) {
  package_dir <- resolve_e8_08_package_dir()

  if (!is.null(package_dir) &&
      requireNamespace("devtools", quietly = TRUE)) {
    load_all <- getExportedValue("devtools", "load_all")
    load_all(package_dir, export_all = FALSE, quiet = TRUE)
  }
}

make_staggered_panel <- function(n_units = 6L, n_periods = 8L, seed = 42L) {
  set.seed(seed)

  panel <- expand.grid(
    id = seq_len(n_units),
    time = seq_len(n_periods),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]

  cohort_map <- c(5, 5, 7, 7, 0, Inf)
  panel$cohort <- cohort_map[panel$id]
  panel$x1 <- rnorm(nrow(panel))
  panel$y <- rnorm(nrow(panel)) +
    ifelse(
      is.finite(panel$cohort) & panel$cohort > 0 & panel$time >= panel$cohort,
      2,
      0
    )

  rownames(panel) <- NULL
  panel
}

make_common_timing_panel <- function(n_units = 10L, seed = 123L) {
  set.seed(seed)

  panel <- expand.grid(
    id = seq_len(n_units),
    time = seq_len(6L),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]

  panel$d <- as.integer(panel$id <= 5L)
  panel$post <- as.integer(panel$time >= 5L)
  panel$x1 <- rnorm(nrow(panel))
  panel$y <- rnorm(nrow(panel)) + ifelse(panel$d == 1L & panel$post == 1L, 1.5, 0)

  rownames(panel) <- NULL
  panel
}

make_unbalanced_panel <- function(seed = 456L) {
  panel <- make_staggered_panel(seed = seed)

  set.seed(seed)
  keep <- sample(
    c(TRUE, FALSE),
    size = nrow(panel),
    replace = TRUE,
    prob = c(0.8, 0.2)
  )

  unbalanced <- panel[keep, , drop = FALSE]
  rownames(unbalanced) <- NULL
  unbalanced
}

make_high_missing_panel <- function(seed = 789L) {
  panel <- make_staggered_panel(seed = seed)

  set.seed(seed)
  missing_index <- sample(
    seq_len(nrow(panel)),
    size = floor(nrow(panel) * 0.35)
  )
  panel$y[missing_index] <- NA_real_

  rownames(panel) <- NULL
  panel
}

make_single_cohort_panel <- function(seed = 101L) {
  set.seed(seed)

  panel <- expand.grid(
    id = seq_len(6L),
    time = seq_len(8L),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]

  panel$cohort <- ifelse(panel$id <= 3L, 5L, 0L)
  panel$y <- rnorm(nrow(panel)) + ifelse(panel$cohort == 5L & panel$time >= 5L, 2, 0)

  rownames(panel) <- NULL
  panel
}

make_e8_08_all_na_outcome_panel <- function() {
  panel <- expand.grid(
    id = seq_len(4L),
    time = seq_len(6L),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]
  panel$cohort <- ifelse(panel$id <= 2L, 5L, 0L)
  panel$y <- NA_real_

  rownames(panel) <- NULL
  panel
}

make_e8_08_inf_nan_outcome_panel <- function() {
  panel <- expand.grid(
    id = seq_len(4L),
    time = seq_len(6L),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]
  panel$cohort <- ifelse(panel$id <= 2L, 5L, 0L)
  panel$y <- seq_len(nrow(panel)) / 10
  panel$y[c(2L, 8L)] <- NaN
  panel$y[c(5L, 9L)] <- Inf
  panel$y[c(6L, 10L)] <- -Inf

  rownames(panel) <- NULL
  panel
}

make_e8_08_single_period_panel <- function() {
  panel <- data.frame(
    id = seq_len(6L),
    time = rep(1L, 6L),
    cohort = c(2L, 2L, 2L, 0L, 0L, 0L),
    y = c(1, 2, 3, 1.5, 2.5, 3.5)
  )

  rownames(panel) <- NULL
  panel
}

make_e8_08_fatal_guard_panel <- function() {
  panel <- expand.grid(
    id = seq_len(6L),
    time = seq_len(4L),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]
  panel$gvar <- c(
    rep(2L, 4L),
    rep(2L, 4L),
    rep(3L, 4L),
    rep(3L, 4L),
    rep(0L, 4L),
    rep(0L, 4L)
  )
  panel$y <- c(
    1.0, 1.2, 1.4, 1.6,
    0.8, 1.0, 1.2, 1.4,
    1.5, 1.7, 1.9, 2.1,
    1.3, 1.5, 1.7, 1.9,
    0.9, 1.1, 1.3, 1.5,
    1.0, 1.2, 1.4, 1.6
  ) + ifelse(panel$gvar > 0 & panel$time >= panel$gvar, 2.0, 0.0)

  rownames(panel) <- NULL
  panel
}

make_seasonal_panel <- function(seed = 202L) {
  set.seed(seed)

  panel <- expand.grid(
    id = seq_len(8L),
    time = seq_len(24L),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]

  panel$cohort <- ifelse(panel$id <= 4L, 13L, 0L)
  panel$y <- rnorm(nrow(panel)) + 2 * sin(2 * pi * panel$time / 4)

  rownames(panel) <- NULL
  panel
}

make_heterogeneous_trend_panel <- function(seed = 303L) {
  set.seed(seed)

  panel <- expand.grid(
    id = seq_len(8L),
    time = seq_len(10L),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]

  panel$cohort <- ifelse(panel$id <= 3L, 6L, ifelse(panel$id <= 6L, 8L, 0L))
  slope <- ifelse(panel$cohort == 6L, 0.5, ifelse(panel$cohort == 8L, 1.5, 0))
  panel$y <- rnorm(nrow(panel)) + slope * panel$time

  rownames(panel) <- NULL
  panel
}

make_clustering_panel <- function(seed = 404L) {
  set.seed(seed)

  panel <- expand.grid(
    id = seq_len(30L),
    time = seq_len(6L),
    KEEP.OUT.ATTRS = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]

  panel$d <- as.integer(panel$id <= 15L)
  panel$post <- as.integer(panel$time >= 4L)
  panel$region <- rep(seq_len(5L), each = 6L)[panel$id]
  panel$state <- rep(seq_len(10L), each = 3L)[panel$id]
  panel$y <- rnorm(nrow(panel)) + ifelse(panel$d == 1L & panel$post == 1L, 1, 0)

  rownames(panel) <- NULL
  panel
}

make_e8_08_short_common_timing_panel <- function() {
  panel <- data.frame(
    firm = rep(1:4, each = 3),
    year = rep(1:3, times = 4),
    y = c(1, 2, 3, 2, 3, 4, 1, 1, 2, 2, 2, 3),
    stringsAsFactors = FALSE
  )

  rownames(panel) <- NULL
  panel
}

make_e8_08_selection_observed_y_panel <- function() {
  data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    gvar = c(2, 2, 0, 0),
    y = c(NA_real_, 5, 3, 4)
  )
}

make_e8_08_selection_attrition_balance_panel <- function() {
  data.frame(
    id = rep(1:4, each = 4),
    time = rep(1:4, times = 4),
    gvar = rep(c(3, 4, 99, 2), each = 4),
    y = c(
      NA_real_, 10, 11, NA_real_,
      5, NA_real_, NA_real_, NA_real_,
      7, 8, 9, 10,
      NA_real_, 9, 10, 11
    )
  )
}

make_e8_08_selection_missing_pattern_cases <- function() {
  list(
    balanced = data.frame(
      id = rep(1:2, each = 3),
      time = rep(1:3, times = 2),
      y = c(1, 2, 3, 4, 5, 6)
    ),
    mcar = data.frame(
      id = rep(1:4, each = 3),
      time = rep(1:3, times = 4),
      y = c(1, 2, 3, 4, NA_real_, 6, 7, 8, 9, 10, 11, 12)
    ),
    mar = data.frame(
      id = rep(1:6, each = 5),
      time = rep(1:5, times = 6),
      gvar = rep(c(4, 4, 5, 0, 0, 0), each = 5),
      x = rep(c(1, 1, 1, 0, 0, 0), each = 5),
      y = c(
        10, 10, NA_real_, 11, 12,
        9, 9, NA_real_, 10, 11,
        8, NA_real_, NA_real_, 9, 10,
        7, 7, 7, 7, 7,
        6, 6, 6, 6, 6,
        5, 5, 5, 5, 5
      )
    ),
    mnar = do.call(
      rbind,
      lapply(seq_len(8L), function(id) {
        data.frame(
          id = id,
          time = 1:4,
          y = if (id <= 4L) {
            c(10, 10, NA_real_, NA_real_)
          } else {
            c(1, 1, 1, 1)
          }
        )
      })
    )
  )
}

make_e8_08_selection_never_treated_filter_panel <- function() {
  data.frame(
    id = rep(1:4, each = 3),
    time = rep(1:3, times = 4),
    gvar = rep(c(3, 0, Inf, NA_real_), each = 3),
    y = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  )
}

make_e8_08_clustering_hierarchy_panel <- function(
    n_regions = 2L,
    states_per_region = 3L,
    units_per_state = 2L,
    periods = 3L
) {
  state_ids <- seq_len(n_regions * states_per_region)
  unit_ids <- seq_len(length(state_ids) * units_per_state)
  cohort_pattern <- rep(c(2005L, 0L, 2007L, 0L, 2005L, 0L), length.out = length(state_ids))

  state_lookup <- data.frame(
    id = unit_ids,
    state = rep(sprintf("s%02d", state_ids), each = units_per_state),
    region = rep(
      sprintf("r%02d", seq_len(n_regions)),
      each = states_per_region * units_per_state
    ),
    cohort = rep(cohort_pattern, each = units_per_state),
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
  panel <- panel[order(panel$id, panel$time), ]

  panel$subcluster <- sprintf("u%02d_t%02d", panel$id, panel$time)
  panel$y <- panel$time + ifelse(panel$cohort > 0L, 0.5, 0)

  rownames(panel) <- NULL
  panel
}

make_e8_08_comprehensive_panel <- function(
    n_units = 40L,
    n_periods = 10L,
    n_treated = 20L,
    treatment_start = 2006L,
    seed = 812L
) {
  set.seed(seed)

  years <- seq.int(treatment_start - (n_periods - 5L), treatment_start + 4L)
  panel <- data.frame(
    id = rep(seq_len(n_units), each = n_periods),
    year = rep(years, times = n_units)
  )

  panel$d <- ifelse(panel$id <= n_treated, 1L, 0L)
  panel$post <- ifelse(panel$year >= treatment_start, 1L, 0L)

  unit_fe <- rnorm(n_units, sd = 0.6)
  x1_unit <- rnorm(n_units)
  x2_unit <- runif(n_units, min = -0.5, max = 0.5)
  time_index <- panel$year - min(panel$year)

  panel$x1 <- x1_unit[panel$id]
  panel$x2 <- x2_unit[panel$id]
  panel$y <- unit_fe[panel$id] +
    0.15 * time_index +
    0.4 * panel$x1 -
    0.25 * panel$x2 +
    1.5 * panel$d * panel$post +
    rnorm(nrow(panel), sd = 0.15)

  rownames(panel) <- NULL
  panel
}

make_e8_08_parallel_trends_common_panel <- function(delta = 0) {
  panel <- expand.grid(
    firm = 1:12,
    year = 1:8,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  panel <- panel[order(panel$firm, panel$year), ]
  panel$eps <- ((panel$firm * 7 + panel$year * 3) %% 11 - 5) / 10
  panel$y <- with(
    panel,
    2 + firm / 5 + 0.4 * year + eps +
      ifelse(firm <= 6 & year < 4, delta * (year - 1), 0) +
      ifelse(firm <= 6 & year >= 4, 1.0, 0)
  )

  rownames(panel) <- NULL
  panel
}

count_e8_08_comprehensive_plot_panels <- function(plot_object) {
  if ("patchwork" %in% class(plot_object)) {
    return(length(plot_object$patches$plots) + 1L)
  }

  if (is.list(plot_object) && !inherits(plot_object, "ggplot")) {
    return(length(plot_object))
  }

  if (inherits(plot_object, c("gtable", "gTree"))) {
    return(2L)
  }

  1L
}

make_e8_08_pre_period_spec <- function(
    spec_id,
    n_pre_periods,
    excluded_periods = 0L,
    att = NA_real_,
    se = 0.10,
    pvalue = 0.01,
    converged = TRUE,
    spec_warnings = character(0)
) {
  if (!converged) {
    att <- NA_real_
    se <- NA_real_
    pvalue <- NA_real_
  }

  t_stat <- if (isTRUE(converged) && !is.na(att) && !is.na(se) && se != 0) {
    att / se
  } else {
    NA_real_
  }

  ci_lower <- if (isTRUE(converged) && !is.na(att) && !is.na(se)) {
    att - 1.96 * se
  } else {
    NA_real_
  }

  ci_upper <- if (isTRUE(converged) && !is.na(att) && !is.na(se)) {
    att + 1.96 * se
  } else {
    NA_real_
  }

  list(
    specification_id = as.integer(spec_id),
    n_pre_periods = as.integer(n_pre_periods),
    start_period = as.integer(0L),
    end_period = as.integer(0L),
    excluded_periods = as.integer(excluded_periods),
    att = att,
    se = se,
    t_stat = t_stat,
    pvalue = pvalue,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    n_treated = if (isTRUE(converged)) 5L else 0L,
    n_control = if (isTRUE(converged)) 5L else 0L,
    df = if (isTRUE(converged)) 8L else 0L,
    converged = converged,
    spec_warnings = spec_warnings
  )
}

resolve_e8_08_parity_path <- function(filename) {
  candidates <- c(
    file.path(
      "/Users/cxy/Desktop/lwdid_r",
      "_automation", "test-artifacts", "parity", filename
    ),
    testthat::test_path(
      "..", "..", "..",
      "_automation", "test-artifacts", "parity", filename
    )
  )
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  existing <- candidates[file.exists(candidates)]

  if (length(existing) > 0L) {
    return(existing[[1L]])
  }

  candidates[[1L]]
}

read_e8_08_parity_json <- function(filename) {
  path <- resolve_e8_08_parity_path(filename)
  if (!file.exists(path)) {
    testthat::skip(paste("parity json is unavailable in this test environment:", path))
  }

  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

read_e8_08_parity_csv <- function(filename) {
  path <- resolve_e8_08_parity_path(filename)
  if (!file.exists(path)) {
    testthat::skip(paste("parity fixture is unavailable in this test environment:", path))
  }

  utils::read.csv(path, stringsAsFactors = FALSE)
}

read_e8_08_required_parity_json <- function(filename) {
  skip_if_not_installed("jsonlite")

  path <- resolve_e8_08_parity_path(filename)
  expect_true(file.exists(path))
  if (!file.exists(path)) {
    return(NULL)
  }

  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

compute_e8_08_smoking_summary <- function() {
  data("smoking", package = "lwdid", envir = environment())

  warnings_seen <- character(0)
  result <- withCallingHandlers(
    lwdid_sensitivity(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      type = "all",
      verbose = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    warnings = sort(unique(warnings_seen)),
    robustness_level = as.character(result$pre_period_result$robustness_level),
    sensitivity_ratio = unname(result$pre_period_result$sensitivity_ratio),
    demean_att = unname(result$transformation_comparison$demean_att)
  )
}

compute_e8_08_readme_contract <- function() {
  data("smoking", package = "lwdid", envir = environment())

  warnings_seen <- character(0)
  run_with_warnings <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }

  pre_period <- run_with_warnings(
    lwdid_sensitivity_pre_period(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      verbose = FALSE
    )
  )
  anticipation <- run_with_warnings(
    sensitivity_no_anticipation(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      verbose = FALSE
    )
  )
  comprehensive <- run_with_warnings(
    lwdid_sensitivity(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      type = "all",
      verbose = FALSE
    )
  )

  list(
    warnings = sort(unique(enc2utf8(warnings_seen))),
    pre_period = list(
      sensitivity_ratio = unname(as.numeric(pre_period$sensitivity_ratio)),
      robustness_level = as.character(pre_period$robustness_level),
      is_robust = isTRUE(pre_period$is_robust)
    ),
    anticipation = list(
      anticipation_detected = isTRUE(anticipation$anticipation_detected),
      recommended_exclusion = unname(as.integer(anticipation$recommended_exclusion)),
      detection_method = as.character(anticipation$detection_method)
    ),
    comprehensive = list(
      overall_assessment = enc2utf8(comprehensive$overall_assessment),
      recommendations = enc2utf8(comprehensive$recommendations),
      transformation = list(
        demean_att = unname(as.numeric(comprehensive$transformation_comparison$demean_att)),
        demean_se = unname(as.numeric(comprehensive$transformation_comparison$demean_se)),
        detrend_att = unname(as.numeric(comprehensive$transformation_comparison$detrend_att)),
        detrend_se = unname(as.numeric(comprehensive$transformation_comparison$detrend_se)),
        difference = unname(as.numeric(comprehensive$transformation_comparison$difference)),
        rel_diff = unname(as.numeric(comprehensive$transformation_comparison$rel_diff))
      ),
      estimator_is_null = is.null(comprehensive$estimator_comparison),
      print_output = enc2utf8(capture.output(print(comprehensive))),
      summary_output = enc2utf8(capture.output(summary(comprehensive)))
    )
  )
}

compute_e8_08_clustering_hierarchical_summary <- function() {
  fixture_path <- resolve_e8_08_parity_path(
    "e8_04_clustering_layer3_hierarchical_fixture.csv"
  )
  fixture <- utils::read.csv(fixture_path, stringsAsFactors = FALSE)

  diagnosis <- diagnose_clustering(
    fixture,
    ivar = "idcode",
    potential_cluster_vars = c("state", "region"),
    gvar = "first_treat",
    d = NULL,
    verbose = FALSE
  )
  recommendation <- recommend_clustering(
    fixture,
    ivar = "idcode",
    tvar = "year",
    potential_cluster_vars = c("state", "region"),
    gvar = "first_treat",
    d = NULL,
    verbose = FALSE
  )
  consistency <- check_clustering_consistency(
    fixture,
    ivar = "idcode",
    cluster_var = "state",
    gvar = "first_treat",
    d = NULL,
    verbose = FALSE
  )

  list(
    diagnosis = list(
      recommended_cluster_var = diagnosis$recommended_cluster_var,
      treatment_variation_level = diagnosis$treatment_variation_level,
      recommendation_reason = diagnosis$recommendation_reason,
      warnings = unname(as.character(diagnosis$warnings))
    ),
    recommendation = list(
      recommended_var = recommendation$recommended_var,
      n_clusters = as.integer(recommendation$n_clusters),
      n_treated_clusters = as.integer(recommendation$n_treated_clusters),
      n_control_clusters = as.integer(recommendation$n_control_clusters),
      confidence = as.numeric(recommendation$confidence),
      reasons = unname(as.character(recommendation$reasons)),
      alternatives = lapply(
        recommendation$alternatives,
        function(option) {
          list(
            var = option$var,
            n_clusters = as.integer(option$n_clusters),
            reliability_score = as.numeric(option$reliability_score),
            reason = option$reason
          )
        }
      ),
      warnings = unname(as.character(recommendation$warnings)),
      use_wild_bootstrap = isTRUE(recommendation$use_wild_bootstrap),
      wild_bootstrap_reason = recommendation$wild_bootstrap_reason
    ),
    consistency = list(
      is_consistent = isTRUE(consistency$is_consistent),
      treatment_variation_level = consistency$treatment_variation_level,
      cluster_level = consistency$cluster_level,
      n_clusters = as.integer(consistency$n_clusters),
      n_inconsistent = as.integer(consistency$n_inconsistent),
      pct_inconsistent = as.numeric(consistency$pct_inconsistent),
      recommendation = consistency$recommendation,
      details = consistency$details
    )
  )
}

compute_e8_08_fatal_guard_summary <- function() {
  panel <- make_e8_08_fatal_guard_panel()

  strict_result <- suppressWarnings(lwdid(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    control_group = "not_yet_treated",
    aggregate = "none",
    verbose = "quiet"
  ))

  strict_rows <- strict_result$cohort_time_effects
  pre_boundary_row <- strict_rows[
    strict_rows$cohort == 2L & strict_rows$period == 2L,
    ,
    drop = FALSE
  ]
  boundary_row <- strict_rows[
    strict_rows$cohort == 2L & strict_rows$period == 3L,
    ,
    drop = FALSE
  ]

  pre_boundary_ids <- panel[
    panel$time == 2L & (panel$gvar > 2L | panel$gvar == 0L),
    "id"
  ]
  boundary_ids <- panel[
    panel$time == 3L & (panel$gvar > 3L | panel$gvar == 0L),
    "id"
  ]

  warning_classes <- character(0)
  warning_messages <- character(0)
  data("castle", package = "lwdid", envir = environment())
  aggregation_result <- withCallingHandlers(
    lwdid(
      castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      aggregate = "overall",
      control_group = "not_yet_treated",
      verbose = "quiet"
    ),
    warning = function(w) {
      warning_classes <<- c(warning_classes, class(w)[[1L]])
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  cohort_n_control <- unname(as.integer(vapply(
    aggregation_result$cohort_effects,
    function(x) x$n_control,
    integer(1)
  )))
  n_never_treated <- as.integer(aggregation_result$n_never_treated)

  list(
    fatal_001 = list(
      focal_cohort = 2L,
      pre_boundary_period = 2L,
      boundary_period = 3L,
      pre_boundary_control_ids = unname(as.integer(pre_boundary_ids)),
      boundary_control_ids = unname(as.integer(boundary_ids)),
      pre_boundary_n_control = as.integer(pre_boundary_row$n_control[[1L]]),
      boundary_n_control = as.integer(boundary_row$n_control[[1L]]),
      pre_boundary_matches_expected =
        identical(as.integer(pre_boundary_row$n_control[[1L]]), length(pre_boundary_ids)),
      boundary_matches_expected =
        identical(as.integer(boundary_row$n_control[[1L]]), length(boundary_ids))
    ),
    fatal_004 = list(
      requested_control_group = as.character(aggregation_result$control_group),
      used_control_group = as.character(aggregation_result$control_group_used),
      n_never_treated = n_never_treated,
      n_control_summary = as.integer(aggregation_result$n_control),
      overall_n_control = as.integer(aggregation_result$overall_effect$n_control),
      cohort_n_control = cohort_n_control,
      warning_classes = unname(unique(warning_classes)),
      warning_messages = unname(warning_messages),
      switch_warning_detected = any(warning_classes == "lwdid_control_group_switch"),
      nt_only_confirmed = identical(
        c(
          as.integer(aggregation_result$n_control),
          as.integer(aggregation_result$overall_effect$n_control),
          cohort_n_control
        ),
        rep.int(n_never_treated, length(cohort_n_control) + 2L)
      )
    )
  )
}

compute_e8_08_return_diagnostics_surface_summary <- function() {
  data("castle", package = "lwdid", envir = environment())
  data("smoking", package = "lwdid", envir = environment())

  castle_diagnosed <- suppressWarnings(lwdid(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_var = "state_name",
    aggregate = "cohort",
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  smoking_messages <- character(0)
  smoking_output <- capture.output(
    smoking_diagnosed <- withCallingHandlers(
      suppressWarnings(lwdid(
        smoking,
        y = "lcigsale",
        ivar = "state",
        tvar = "year",
        d = "d",
        post = "post",
        return_diagnostics = TRUE,
        verbose = "quiet"
      )),
      message = function(m) {
        smoking_messages <<- c(smoking_messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )

  list(
    castle = list(
      diagnostics_names = names(castle_diagnosed[["diagnostics"]]),
      selection_class = class(castle_diagnosed[["diagnostics"]][["selection"]])[[1L]],
      trends_class = class(castle_diagnosed[["diagnostics"]][["trends"]])[[1L]],
      parallel_trends_class =
        class(castle_diagnosed[["diagnostics"]][["parallel_trends"]])[[1L]],
      clustering_class = class(castle_diagnosed[["diagnostics"]][["clustering"]])[[1L]],
      controls_tier_is_null =
        is.null(castle_diagnosed[["diagnostics"]][["controls_tier"]])
    ),
    smoking = list(
      diagnostics_names = names(smoking_diagnosed[["diagnostics"]]),
      controls_tier =
        as.character(smoking_diagnosed[["diagnostics"]][["controls_tier"]]),
      selection_class = class(smoking_diagnosed[["diagnostics"]][["selection"]])[[1L]],
      trends_class = class(smoking_diagnosed[["diagnostics"]][["trends"]])[[1L]],
      parallel_trends_class =
        class(smoking_diagnosed[["diagnostics"]][["parallel_trends"]])[[1L]],
      clustering_is_null =
        is.null(smoking_diagnosed[["diagnostics"]][["clustering"]]),
      captured_messages = unname(enc2utf8(smoking_messages)),
      captured_output = unname(enc2utf8(smoking_output))
    )
  )
}

compute_e8_08_diagnostics_invariance_summary <- function() {
  data("castle", package = "lwdid", envir = environment())
  data("smoking", package = "lwdid", envir = environment())

  castle_base <- suppressWarnings(lwdid(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_var = "state_name",
    aggregate = "cohort",
    verbose = "quiet"
  ))
  castle_diagnosed <- suppressWarnings(lwdid(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_var = "state_name",
    aggregate = "cohort",
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  smoking_base <- suppressWarnings(lwdid(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    post = "post",
    verbose = "quiet"
  ))

  smoking_messages <- character(0)
  smoking_output <- capture.output(
    smoking_diagnosed <- withCallingHandlers(
      suppressWarnings(lwdid(
        smoking,
        y = "lcigsale",
        ivar = "state",
        tvar = "year",
        d = "d",
        post = "post",
        return_diagnostics = TRUE,
        verbose = "quiet"
      )),
      message = function(m) {
        smoking_messages <<- c(smoking_messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )

  list(
    castle = list(
      base = list(
        att = unname(as.numeric(castle_base$att)),
        se_att = unname(as.numeric(castle_base$se_att)),
        df_inference = unname(as.integer(castle_base$df_inference))
      ),
      diagnosed = list(
        att = unname(as.numeric(castle_diagnosed$att)),
        se_att = unname(as.numeric(castle_diagnosed$se_att)),
        df_inference = unname(as.integer(castle_diagnosed$df_inference)),
        diagnostics = list(
          selection_class = class(get_diagnostics(castle_diagnosed, "selection"))[[1L]],
          trends_class = class(get_diagnostics(castle_diagnosed, "trends"))[[1L]],
          clustering_class = class(get_diagnostics(castle_diagnosed, "clustering"))[[1L]],
          controls_tier = if (is.null(castle_diagnosed$diagnostics$controls_tier)) {
            NA_character_
          } else {
            as.character(castle_diagnosed$diagnostics$controls_tier)
          }
        )
      )
    ),
    smoking = list(
      base = list(
        att = unname(as.numeric(smoking_base$att)),
        se_att = unname(as.numeric(smoking_base$se_att)),
        df_inference = unname(as.integer(smoking_base$df_inference))
      ),
      diagnosed = list(
        att = unname(as.numeric(smoking_diagnosed$att)),
        se_att = unname(as.numeric(smoking_diagnosed$se_att)),
        df_inference = unname(as.integer(smoking_diagnosed$df_inference)),
        diagnostics = list(
          selection_class = class(get_diagnostics(smoking_diagnosed, "selection"))[[1L]],
          trends_class = class(get_diagnostics(smoking_diagnosed, "trends"))[[1L]],
          clustering_is_null = is.null(get_diagnostics(smoking_diagnosed, "clustering")),
          controls_tier = if (is.null(smoking_diagnosed$diagnostics$controls_tier)) {
            NA_character_
          } else {
            as.character(smoking_diagnosed$diagnostics$controls_tier)
          }
        ),
        captured_messages = unname(enc2utf8(smoking_messages)),
        captured_output = unname(enc2utf8(smoking_output))
      )
    )
  )
}

compute_e8_08_return_diagnostics_surface_summary <- function() {
  data("castle", package = "lwdid", envir = environment())
  data("smoking", package = "lwdid", envir = environment())

  castle_diagnosed <- suppressWarnings(lwdid(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_var = "state_name",
    aggregate = "cohort",
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  smoking_messages <- character(0)
  smoking_output <- capture.output(
    smoking_diagnosed <- withCallingHandlers(
      suppressWarnings(lwdid(
        smoking,
        y = "lcigsale",
        ivar = "state",
        tvar = "year",
        d = "d",
        post = "post",
        return_diagnostics = TRUE,
        verbose = "quiet"
      )),
      message = function(m) {
        smoking_messages <<- c(smoking_messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )

  list(
    castle = list(
      diagnostics_names = names(castle_diagnosed$diagnostics),
      selection_class = class(castle_diagnosed$diagnostics$selection)[[1L]],
      trends_class = class(castle_diagnosed$diagnostics$trends)[[1L]],
      parallel_trends_class = class(castle_diagnosed$diagnostics$parallel_trends)[[1L]],
      clustering_class = class(castle_diagnosed$diagnostics$clustering)[[1L]],
      controls_tier_is_null = is.null(castle_diagnosed$diagnostics$controls_tier)
    ),
    smoking = list(
      diagnostics_names = names(smoking_diagnosed$diagnostics),
      controls_tier = if (is.null(smoking_diagnosed$diagnostics$controls_tier)) {
        NA_character_
      } else {
        as.character(smoking_diagnosed$diagnostics$controls_tier)
      },
      selection_class = class(smoking_diagnosed$diagnostics$selection)[[1L]],
      trends_class = class(smoking_diagnosed$diagnostics$trends)[[1L]],
      parallel_trends_class = class(smoking_diagnosed$diagnostics$parallel_trends)[[1L]],
      clustering_is_null = is.null(smoking_diagnosed$diagnostics$clustering),
      captured_messages = unname(enc2utf8(smoking_messages)),
      captured_output = unname(enc2utf8(smoking_output))
    )
  )
}

compute_e8_08_get_diagnostics_summary <- function() {
  data("castle", package = "lwdid", envir = environment())
  data("smoking", package = "lwdid", envir = environment())

  castle_diagnosed <- suppressWarnings(lwdid(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_var = "state_name",
    aggregate = "cohort",
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  smoking_messages <- character(0)
  smoking_output <- capture.output(
    smoking_diagnosed <- withCallingHandlers(
      suppressWarnings(lwdid(
        smoking,
        y = "lcigsale",
        ivar = "state",
        tvar = "year",
        d = "d",
        post = "post",
        return_diagnostics = TRUE,
        verbose = "quiet"
      )),
      message = function(m) {
        smoking_messages <<- c(smoking_messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )

  castle_all <- get_diagnostics(castle_diagnosed, "all")
  smoking_all <- get_diagnostics(smoking_diagnosed, "all")

  list(
    castle = list(
      all_names = names(castle_all),
      selection_class = class(get_diagnostics(castle_diagnosed, "selection"))[[1L]],
      trends_class = class(get_diagnostics(castle_diagnosed, "trends"))[[1L]],
      clustering_class = class(get_diagnostics(castle_diagnosed, "clustering"))[[1L]],
      sensitivity_is_null = is.null(get_diagnostics(castle_diagnosed, "sensitivity"))
    ),
    smoking = list(
      all_names = names(smoking_all),
      selection_class = class(get_diagnostics(smoking_diagnosed, "selection"))[[1L]],
      trends_class = class(get_diagnostics(smoking_diagnosed, "trends"))[[1L]],
      clustering_is_null = is.null(get_diagnostics(smoking_diagnosed, "clustering")),
      sensitivity_is_null = is.null(get_diagnostics(smoking_diagnosed, "sensitivity")),
      captured_messages = unname(enc2utf8(smoking_messages)),
      captured_output = unname(enc2utf8(smoking_output))
    )
  )
}

compute_e8_08_default_get_diagnostics_summary <- function() {
  data("smoking", package = "lwdid", envir = environment())

  result <- suppressWarnings(lwdid(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    post = "post",
    verbose = "quiet"
  ))

  captured_messages <- character(0)
  returned <- withCallingHandlers(
    get_diagnostics(result),
    message = function(m) {
      captured_messages <<- c(captured_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  diagnostics <- result[["diagnostics"]]

  list(
    diagnostics_names = names(diagnostics),
    controls_tier = as.character(diagnostics[["controls_tier"]]),
    returned_is_null = is.null(returned),
    captured_messages = unname(enc2utf8(captured_messages))
  )
}

compute_e8_08_diagnose_suite_summary <- function() {
  data("castle", package = "lwdid", envir = environment())
  data("smoking", package = "lwdid", envir = environment())

  castle_suite <- suppressWarnings(lwdid_diagnose(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_vars = c("state_name"),
    verbose = FALSE
  ))
  smoking_suite <- suppressWarnings(lwdid_diagnose(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    controls = c("lretprice", "lnincome"),
    cluster_vars = NULL,
    verbose = FALSE
  ))

  list(
    castle = list(
      class = class(castle_suite),
      names = names(castle_suite),
      selection_class = class(castle_suite$selection)[[1L]],
      trends_class = class(castle_suite$trends)[[1L]],
      clustering_class = class(castle_suite$clustering)[[1L]],
      recommended_cluster_var = castle_suite$clustering$recommended_cluster_var,
      treatment_variation_level = castle_suite$clustering$treatment_variation_level,
      selection_risk = castle_suite$selection$selection_risk,
      missing_rate_overall = castle_suite$selection$missing_rate_overall,
      recommended_method = castle_suite$trends$recommended_method,
      confidence_level = castle_suite$trends$confidence_level
    ),
    smoking = list(
      class = class(smoking_suite),
      names = names(smoking_suite),
      selection_class = class(smoking_suite$selection)[[1L]],
      trends_class = class(smoking_suite$trends)[[1L]],
      clustering_is_null = is.null(smoking_suite$clustering),
      selection_risk = smoking_suite$selection$selection_risk,
      missing_rate_overall = smoking_suite$selection$missing_rate_overall,
      recommended_method = smoking_suite$trends$recommended_method,
      confidence_level = smoking_suite$trends$confidence_level
    )
  )
}

compute_e8_08_diagnose_failure_isolation_summary <- function() {
  data("castle", package = "lwdid", envir = environment())

  caught_warnings <- character(0)
  suite <- withCallingHandlers(
    lwdid_diagnose(
      castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      controls = c("income", "unemployrt", "poverty"),
      cluster_vars = c("missing_cluster"),
      verbose = FALSE
    ),
    warning = function(w) {
      caught_warnings <<- c(caught_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    class = class(suite),
    names = names(suite),
    selection_class = class(suite$selection)[[1L]],
    trends_class = class(suite$trends)[[1L]],
    clustering_is_null = is.null(suite$clustering),
    selection_risk = suite$selection$selection_risk,
    recommended_method = suite$trends$recommended_method,
    confidence_level = suite$trends$confidence_level,
    warnings_have_missing_cluster =
      any(grepl("missing_cluster", caught_warnings, fixed = TRUE))
  )
}

compute_e8_08_diagnose_print_summary_summary <- function() {
  data("castle", package = "lwdid", envir = environment())
  data("smoking", package = "lwdid", envir = environment())

  castle_suite <- suppressWarnings(lwdid_diagnose(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_vars = c("state_name"),
    verbose = FALSE
  ))
  smoking_suite <- suppressWarnings(lwdid_diagnose(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    controls = c("lretprice", "lnincome"),
    cluster_vars = NULL,
    verbose = FALSE
  ))

  list(
    castle = list(
      print_output = enc2utf8(capture.output(print(castle_suite))),
      summary_output = enc2utf8(capture.output(summary(castle_suite)))
    ),
    smoking = list(
      print_output = enc2utf8(capture.output(print(smoking_suite))),
      summary_output = enc2utf8(capture.output(summary(smoking_suite)))
    )
  )
}

compute_e8_08_never_treated_propagation_summary <- function() {
  data("castle", package = "lwdid", envir = environment())

  late_sid <- unique(castle$sid[castle$gvar == max(castle$gvar, na.rm = TRUE)])
  marker_zero <- castle
  marker_zero$gvar[marker_zero$sid %in% late_sid] <- 0

  marker_custom <- marker_zero
  marker_custom$gvar[marker_custom$gvar == 0] <- 999

  zero_suite <- suppressWarnings(lwdid_diagnose(
    marker_zero,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_vars = NULL,
    verbose = FALSE
  ))
  custom_suite <- suppressWarnings(lwdid_diagnose(
    marker_custom,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_vars = NULL,
    never_treated_values = 999,
    verbose = FALSE
  ))

  summarize_suite <- function(suite) {
    list(
      class = class(suite)[[1L]],
      selection_class = class(suite$selection)[[1L]],
      trends_class = class(suite$trends)[[1L]],
      clustering_is_null = is.null(suite$clustering),
      selection_risk = suite$selection$selection_risk,
      missing_rate_overall = unname(as.numeric(suite$selection$missing_rate_overall)),
      recommended_method = suite$trends$recommended_method,
      confidence_level = suite$trends$confidence_level
    )
  }

  list(
    zero_marker = summarize_suite(zero_suite),
    custom_marker = summarize_suite(custom_suite)
  )
}

compute_e8_08_diagnostics_plot_summary <- function() {
  data("castle", package = "lwdid", envir = environment())
  data("smoking", package = "lwdid", envir = environment())

  smoking_result <- suppressWarnings(lwdid(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    post = "post",
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))
  castle_result <- suppressWarnings(lwdid(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_var = "state_name",
    aggregate = "cohort",
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  pluck_label <- function(plot_obj, name) {
    label <- plot_obj$labels[[name]]
    if (is.null(label)) {
      return(NA_character_)
    }
    as.character(label)
  }

  summarize_selection_plot <- function(plot_obj) {
    list(
      class = class(plot_obj),
      inherits_ggplot = inherits(plot_obj, "ggplot"),
      inherits_patchwork = inherits(plot_obj, "patchwork"),
      title = pluck_label(plot_obj, "title"),
      x_label = pluck_label(plot_obj, "x"),
      y_label = pluck_label(plot_obj, "y")
    )
  }

  summarize_panel_plot <- function(plot_obj) {
    patch_list <- plot_obj$patches$plots
    if (is.null(patch_list)) {
      patch_list <- list()
    }

    annotation_title <- plot_obj$patches$annotation$title
    if (is.null(annotation_title)) {
      annotation_title <- NA_character_
    } else {
      annotation_title <- as.character(annotation_title)
    }

    patch_titles <- unname(vapply(
      patch_list,
      function(subplot) pluck_label(subplot, "title"),
      character(1)
    ))

    list(
      class = class(plot_obj),
      inherits_ggplot = inherits(plot_obj, "ggplot"),
      inherits_patchwork = inherits(plot_obj, "patchwork"),
      base_title = pluck_label(plot_obj, "title"),
      annotation_title = annotation_title,
      patch_count = length(patch_list),
      patch_titles = patch_titles
    )
  }

  list(
    smoking_selection = summarize_selection_plot(plot(
      smoking_result,
      type = "diagnostics",
      which = "selection"
    )),
    smoking_panel = summarize_panel_plot(plot(
      smoking_result,
      type = "diagnostics"
    )),
    castle_panel = summarize_panel_plot(plot(
      castle_result,
      type = "diagnostics"
    ))
  )
}

compute_e8_08_monte_carlo_summary <- function(fixture_df) {
  split_fixture <- split(fixture_df, fixture_df$replication)

  rows <- lapply(split_fixture, function(rep_df) {
    result <- suppressMessages(withCallingHandlers(
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
        invokeRestart("muffleWarning")
      }
    ))

    demean_att <- result$transformation_comparison$demean_att
    demean_se <- result$transformation_comparison$demean_se
    true_att <- unique(rep_df$true_att)
    ci_lower <- demean_att - 1.96 * demean_se
    ci_upper <- demean_att + 1.96 * demean_se

    data.frame(
      demean_att = demean_att,
      demean_se = demean_se,
      true_att = true_att,
      covered = ci_lower <= true_att && true_att <= ci_upper
    )
  })

  metrics <- do.call(rbind, rows)
  metrics$covered <- as.logical(metrics$covered)

  list(
    mean_att = mean(metrics$demean_att),
    mean_se = mean(metrics$demean_se),
    coverage = mean(metrics$covered)
  )
}

test_that("E8-08.1 scaffold: staggered fixture builder is available and reproducible", {
  first <- make_staggered_panel()
  second <- make_staggered_panel()

  expect_identical(first, second)
  expect_equal(nrow(first), 48L)
  expect_equal(sort(unique(first$cohort)), c(0, 5, 7, Inf))
})

test_that("E8-08.1 scaffold: common-timing fixture builder is available and reproducible", {
  first <- make_common_timing_panel()
  second <- make_common_timing_panel()

  expect_identical(first, second)
  expect_equal(nrow(first), 60L)
  expect_equal(sum(first$post == 1L), 20L)
})

test_that("E8-08.1 scaffold: missingness-oriented fixture builders create targeted panels", {
  unbalanced <- make_unbalanced_panel()
  high_missing <- make_high_missing_panel()

  expect_identical(unbalanced, make_unbalanced_panel())
  expect_identical(high_missing, make_high_missing_panel())
  expect_lt(nrow(unbalanced), nrow(make_staggered_panel()))
  expect_gt(mean(is.na(high_missing$y)), 0.30)
})

test_that("E8-08.1 scaffold: trend-oriented fixture builders expose single-cohort and heterogeneity structure", {
  single_cohort <- make_single_cohort_panel()
  seasonal <- make_seasonal_panel()
  heterogeneous <- make_heterogeneous_trend_panel()

  expect_identical(single_cohort, make_single_cohort_panel())
  expect_identical(seasonal, make_seasonal_panel())
  expect_identical(heterogeneous, make_heterogeneous_trend_panel())
  expect_equal(sort(unique(single_cohort$cohort)), c(0, 5))
  expect_equal(max(seasonal$time), 24L)

  cohort_means <- aggregate(y ~ cohort + time, data = heterogeneous, FUN = mean)
  slope_6 <- with(subset(cohort_means, cohort == 6), y[time == 10] - y[time == 1])
  slope_8 <- with(subset(cohort_means, cohort == 8), y[time == 10] - y[time == 1])
  expect_gt(slope_8, slope_6)
})

test_that("E8-08.1 scaffold: clustering fixture builder carries candidate cluster variables", {
  clustering <- make_clustering_panel()

  expect_identical(clustering, make_clustering_panel())
  expect_equal(length(unique(clustering$id)), 30L)
  expect_equal(length(unique(clustering$region)), 5L)
  expect_equal(length(unique(clustering$state)), 10L)
  expect_true(all(c("d", "post", "y") %in% names(clustering)))
})

test_that("E8-08 A1 helper: sensitivity ratio reproduces Python-specified arithmetic and edge guards", {
  oracle <- read_e8_08_required_parity_json(
    "20260327-qa-parity-e8-08-task9-numeric-oracle.json"
  )
  if (is.null(oracle)) {
    return(invisible(NULL))
  }

  expect_equal(
    lwdid:::.compute_sensitivity_ratio(c(1.0, 1.5, 2.0), 1.5),
    (2.0 - 1.0) / abs(1.5),
    tolerance = 1e-10
  )
  expect_equal(
    lwdid:::.compute_sensitivity_ratio(c(-3.0, -1.0), -2.0),
    1.0,
    tolerance = 1e-10
  )
  expect_equal(
    lwdid:::.compute_sensitivity_ratio(numeric(0), 1.0),
    0,
    tolerance = 1e-10
  )
  expect_equal(lwdid:::.compute_sensitivity_ratio(c(0.1, 0.2), 1e-15), Inf)

  sr_case <- oracle$cases$sensitivity_ratio_formula
  current_results <- vapply(
    sr_case$manual_cases,
    function(case) {
      lwdid:::.compute_sensitivity_ratio(
        unlist(case$atts, use.names = FALSE),
        as.numeric(case$baseline_att)
      )
    },
    numeric(1)
  )
  expected_results <- vapply(
    sr_case$manual_cases,
    function(case) as.numeric(case$manual_expected),
    numeric(1)
  )

  expect_equal(
    current_results[["positive_range"]],
    expected_results[["positive_range"]],
    tolerance = 1e-12
  )
  expect_equal(
    current_results[["negative_range"]],
    expected_results[["negative_range"]],
    tolerance = 1e-12
  )
  expect_equal(
    current_results[["empty_att_vector"]],
    expected_results[["empty_att_vector"]],
    tolerance = 1e-12
  )
  expect_true(is.infinite(current_results[["near_zero_baseline"]]))
  expect_identical(sr_case$comparison$status, "matched")
})

test_that("E8-08 A2 helper: sensitivity ratio returns zero for flat ranges including zero baseline", {
  expect_equal(
    lwdid:::.compute_sensitivity_ratio(c(1.0, 1.0, 1.0), 1.0),
    0,
    tolerance = 1e-10
  )
  expect_equal(
    lwdid:::.compute_sensitivity_ratio(c(0, 0), 0),
    0,
    tolerance = 1e-10
  )
})

test_that("E8-08 A3 helper: robustness bins honor exact threshold cutoffs", {
  expect_identical(lwdid:::.determine_robustness_level(0.0999), "highly_robust")
  expect_identical(lwdid:::.determine_robustness_level(0.10), "moderately_robust")
  expect_identical(lwdid:::.determine_robustness_level(0.2499), "moderately_robust")
  expect_identical(lwdid:::.determine_robustness_level(0.25), "sensitive")
  expect_identical(lwdid:::.determine_robustness_level(0.4999), "sensitive")
  expect_identical(lwdid:::.determine_robustness_level(0.50), "highly_sensitive")
})

test_that("E8-08 A4 helper: pre-period auto-detection respects rolling-specific minimums", {
  staggered <- data.frame(
    id = rep(1:6, each = 8),
    time = rep(1:8, times = 6),
    cohort = rep(c(5, 5, 7, 7, 0, Inf), each = 8),
    y = seq_len(48)
  )
  common <- make_common_timing_panel()

  expect_equal(
    lwdid:::.auto_detect_pre_period_range(
      staggered, "id", "time", "cohort", NULL, NULL, "demean"
    ),
    c(1L, 4L)
  )
  expect_equal(
    lwdid:::.auto_detect_pre_period_range(
      staggered, "id", "time", "cohort", NULL, NULL, "detrend"
    ),
    c(2L, 4L)
  )
  expect_equal(
    lwdid:::.auto_detect_pre_period_range(
      common, "id", "time", NULL, "d", "post", "detrendq"
    ),
    c(2L, 4L)
  )
})

test_that("E8-08 A5 api: pre-period sensitivity preserves failed specifications in bookkeeping", {
  panel <- make_common_timing_panel()

  testthat::local_mocked_bindings(
    .run_single_specification = function(..., n_pre_periods, spec_id,
                                         exclude_periods = 0L) {
      switch(
        as.character(n_pre_periods),
        "2" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 2L,
          excluded_periods = exclude_periods,
          att = 1.10,
          pvalue = 0.02
        ),
        "3" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 3L,
          excluded_periods = exclude_periods,
          converged = FALSE,
          spec_warnings = "forced failure"
        ),
        "4" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 4L,
          excluded_periods = exclude_periods,
          att = 1.25,
          pvalue = 0.03
        ),
        stop("unexpected n_pre_periods")
      )
    },
    .package = "lwdid"
  )

  result <- lwdid::lwdid_sensitivity_pre_period(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    d = "d",
    post = "post",
    rolling = "demean",
    pre_period_range = c(2L, 4L),
    verbose = FALSE
  )

  expect_equal(result$n_specifications, 3L)
  expect_equal(length(result$specifications), 3L)
  expect_false(result$specifications[[2L]]$converged)
  expect_true(is.na(result$specifications[[2L]]$att))
  expect_identical(result$specifications[[2L]]$spec_warnings, "forced failure")
  expect_equal(result$baseline_spec$n_pre_periods, 4L)
})

test_that("E8-08 A6 api: pre-period sensitivity tracks sign changes across converged specifications", {
  panel <- make_common_timing_panel()

  testthat::local_mocked_bindings(
    .run_single_specification = function(..., n_pre_periods, spec_id,
                                         exclude_periods = 0L) {
      switch(
        as.character(n_pre_periods),
        "2" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 2L,
          excluded_periods = exclude_periods,
          att = 1.10,
          pvalue = 0.02
        ),
        "3" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 3L,
          excluded_periods = exclude_periods,
          att = -0.40,
          pvalue = 0.02
        ),
        "4" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 4L,
          excluded_periods = exclude_periods,
          att = 1.25,
          pvalue = 0.03
        ),
        stop("unexpected n_pre_periods")
      )
    },
    .package = "lwdid"
  )

  result <- lwdid::lwdid_sensitivity_pre_period(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    d = "d",
    post = "post",
    rolling = "demean",
    pre_period_range = c(2L, 4L),
    verbose = FALSE
  )

  expect_false(result$all_same_sign)
  expect_equal(result$n_sign_changes, 1L)
  expect_equal(result$n_significant, 3L)
})

test_that("E8-08 A7 api: monotonic ATT paths surface a monotonic-trend warning", {
  panel <- make_common_timing_panel()

  testthat::local_mocked_bindings(
    .run_single_specification = function(..., n_pre_periods, spec_id,
                                         exclude_periods = 0L) {
      switch(
        as.character(n_pre_periods),
        "2" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 2L,
          excluded_periods = exclude_periods,
          att = 1.00,
          pvalue = 0.02
        ),
        "3" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 3L,
          excluded_periods = exclude_periods,
          att = 1.50,
          pvalue = 0.02
        ),
        "4" = make_e8_08_pre_period_spec(
          spec_id = spec_id,
          n_pre_periods = 4L,
          excluded_periods = exclude_periods,
          att = 2.00,
          pvalue = 0.02
        ),
        stop("unexpected n_pre_periods")
      )
    },
    .package = "lwdid"
  )

  result <- lwdid::lwdid_sensitivity_pre_period(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    d = "d",
    post = "post",
    rolling = "demean",
    pre_period_range = c(2L, 4L),
    verbose = FALSE
  )

  expect_true(any(grepl("monotonic", result$result_warnings, ignore.case = TRUE)))
})

test_that("E8-08 A8 helper: anticipation helper returns none_detected for stable specifications", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.00, converged = TRUE),
    list(excluded_periods = 1L, att = 1.05, converged = TRUE),
    list(excluded_periods = 2L, att = 1.04, converged = TRUE)
  )

  result <- lwdid:::.detect_anticipation_effects(
    estimates = estimates,
    baseline = estimates[[1L]],
    detection_threshold = 0.10
  )

  expect_false(result$detected)
  expect_identical(result$method, "none_detected")
  expect_equal(result$recommended_exclusion, 0L)
})

test_that("E8-08 A9 helper: coefficient-change detection returns actual excluded_periods", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.00, converged = TRUE),
    list(excluded_periods = 2L, att = 1.25, converged = TRUE),
    list(excluded_periods = 4L, att = 1.12, converged = TRUE)
  )

  result <- lwdid:::.detect_anticipation_effects(
    estimates = estimates,
    baseline = estimates[[1L]],
    detection_threshold = 0.10
  )

  expect_true(result$detected)
  expect_identical(result$method, "coefficient_change")
  expect_equal(result$recommended_exclusion, 2L)
})

test_that("E8-08 A10 helper: trend-break detection preserves actual exclusion counts on gapped paths", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.00, converged = TRUE),
    list(excluded_periods = 2L, att = 1.30, converged = TRUE),
    list(excluded_periods = 4L, att = 1.50, converged = TRUE),
    list(excluded_periods = 5L, att = 1.55, converged = TRUE)
  )

  result <- lwdid:::.detect_anticipation_effects(
    estimates = estimates,
    baseline = estimates[[1L]],
    detection_threshold = 0.80
  )

  expect_true(result$detected)
  expect_identical(result$method, "trend_break")
  expect_equal(result$recommended_exclusion, 4L)
})

test_that("E8-08 A11 helper: direction check blocks false positives when exclusions shrink ATT magnitude", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.00, converged = TRUE),
    list(excluded_periods = 1L, att = 0.70, converged = TRUE),
    list(excluded_periods = 2L, att = 0.72, converged = TRUE)
  )

  result <- lwdid:::.detect_anticipation_effects(
    estimates = estimates,
    baseline = estimates[[1L]],
    detection_threshold = 0.10
  )

  expect_false(result$detected)
  expect_identical(result$method, "none_detected")
  expect_equal(result$recommended_exclusion, 0L)
})

test_that("E8-08 A12 helper: anticipation helper reports insufficient_data with fewer than two valid estimates", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.0, converged = TRUE),
    list(excluded_periods = 1L, att = NA_real_, converged = FALSE)
  )

  result <- lwdid:::.detect_anticipation_effects(
    estimates = estimates,
    baseline = estimates[[1L]],
    detection_threshold = 0.10
  )

  expect_false(result$detected)
  expect_identical(result$method, "insufficient_data")
  expect_equal(result$recommended_exclusion, 0L)
})

test_that("E8-08 A13 api: anticipation analysis caps exclusions at the feasible pre-period range", {
  panel <- make_common_timing_panel()

  result <- suppressWarnings(
    lwdid::sensitivity_no_anticipation(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "d",
      post = "post",
      max_anticipation = 10L,
      verbose = FALSE
    )
  )

  expect_equal(result$max_anticipation_tested, 3L)
  expect_equal(result$n_estimates, 4L)
  expect_equal(
    vapply(result$estimates, function(est) est$excluded_periods, integer(1)),
    0:3
  )
  expect_equal(
    vapply(result$estimates, function(est) est$n_pre_periods_used, integer(1)),
    c(4L, 3L, 2L, 1L)
  )
})

test_that("E8-08 A14 api: comprehensive sensitivity returns all four result blocks", {
  panel <- make_e8_08_comprehensive_panel()

  result <- suppressWarnings(
    lwdid::lwdid_sensitivity(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = c("x1", "x2"),
      type = "all",
      max_anticipation = 2L,
      detection_threshold = 0.90,
      verbose = FALSE
    )
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_setequal(
    names(result),
    c(
      "pre_period_result",
      "no_anticipation_result",
      "transformation_comparison",
      "estimator_comparison",
      "overall_assessment",
      "recommendations"
    )
  )
  expect_false(is.null(result$pre_period_result))
  expect_false(is.null(result$no_anticipation_result))
  expect_false(is.null(result$transformation_comparison))
  expect_false(is.null(result$estimator_comparison))
  expect_equal(result$no_anticipation_result$max_anticipation_tested, 2L)
  expect_equal(result$no_anticipation_result$detection_threshold, 0.90)
})

test_that("E8-08 A15 api: transformation comparison exposes demean and detrend summary fields", {
  panel <- make_e8_08_comprehensive_panel()

  result <- suppressWarnings(
    lwdid::lwdid_sensitivity(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = c("x1", "x2"),
      type = "all",
      max_anticipation = 2L,
      detection_threshold = 0.90,
      verbose = FALSE
    )
  )

  comparison <- result$transformation_comparison
  expected_fields <- c(
    "demean_att",
    "demean_se",
    "detrend_att",
    "detrend_se",
    "difference",
    "rel_diff"
  )

  expect_false(is.null(comparison))
  expect_setequal(names(comparison), expected_fields)
  expect_true(all(vapply(
    expected_fields,
    function(field) is.finite(comparison[[field]]),
    logical(1)
  )))
  expect_equal(
    comparison$difference,
    abs(comparison$demean_att - comparison$detrend_att),
    tolerance = 1e-10
  )
  expect_equal(
    comparison$rel_diff,
    comparison$difference / abs(comparison$demean_att),
    tolerance = 1e-10
  )
})

test_that("E8-08 A17 api: comprehensive plot only assembles pre-period and anticipation panels", {
  skip_if_not_installed("ggplot2")

  panel <- make_e8_08_comprehensive_panel()

  result <- suppressWarnings(
    lwdid::lwdid_sensitivity(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = c("x1", "x2"),
      type = "all",
      max_anticipation = 2L,
      detection_threshold = 0.90,
      verbose = FALSE
    )
  )

  plot_object <- plot(result)

  expect_equal(count_e8_08_comprehensive_plot_panels(plot_object), 2L)
  expect_false(is.null(result$transformation_comparison))
  expect_false(is.null(result$estimator_comparison))
})

test_that("E8-08 A18 helper: zero-baseline transformation guard avoids false sensitivity flags", {
  near_zero <- lwdid:::.summarize_transformation_comparison(
    demean_estimate = list(att = 1e-11, se_att = 0.2),
    detrend_estimate = list(att = 2e-11, se_att = 0.3)
  )
  flagged <- lwdid:::.summarize_transformation_comparison(
    demean_estimate = list(att = 1e-9, se_att = 0.2),
    detrend_estimate = list(att = 2e-9, se_att = 0.3)
  )

  expect_equal(near_zero$difference, 1e-11)
  expect_true(is.na(near_zero$rel_diff))
  expect_equal(flagged$difference, 1e-9)
  expect_equal(flagged$rel_diff, 1)

  zero_assessment <- lwdid:::.compute_overall_assessment(
    pre_period_result = NULL,
    no_anticipation_result = NULL,
    transformation_comparison = near_zero,
    estimator_comparison = NULL
  )
  flagged_assessment <- lwdid:::.compute_overall_assessment(
    pre_period_result = NULL,
    no_anticipation_result = NULL,
    transformation_comparison = flagged,
    estimator_comparison = NULL
  )

  expect_identical(
    zero_assessment$overall_assessment,
    "Estimates are robust across multiple robustness checks"
  )
  expect_identical(
    zero_assessment$recommendations,
    "No major robustness issues found"
  )
  expect_identical(
    flagged_assessment$overall_assessment,
    "Caution: transformation sensitivity detected, see recommendations"
  )
  expect_true(any(grepl("detrend", flagged_assessment$recommendations, fixed = TRUE)))
})

test_that("E8-08 A19 api: single specification keeps att_std at zero instead of NA", {
  panel <- make_common_timing_panel()

  result <- suppressWarnings(
    lwdid::lwdid_sensitivity_pre_period(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "d",
      post = "post",
      rolling = "demean",
      pre_period_range = c(4L, 4L),
      verbose = FALSE
    )
  )

  expect_equal(result$n_specifications, 1L)
  expect_equal(result$att_std, 0.0, tolerance = 1e-10)
  expect_false(is.na(result$att_std))
})

test_that("E8-08 A16 api: estimator comparison carries finite ATT and SE per estimator", {
  panel <- make_common_timing_panel()
  unit_controls <- seq_along(unique(panel$id)) / 10
  panel$x1 <- unit_controls[panel$id]

  result <- suppressWarnings(
    lwdid::lwdid_sensitivity(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "d",
      post = "post",
      controls = "x1",
      type = "estimator",
      verbose = FALSE
    )
  )

  estimator_comparison <- result$estimator_comparison
  expected_fields <- c("ra", "ra_se", "ipw", "ipw_se", "ipwra", "ipwra_se")

  expect_false(is.null(estimator_comparison))
  expect_true(all(expected_fields %in% names(estimator_comparison)))
  expect_true(all(vapply(
    expected_fields,
    function(field) is.finite(estimator_comparison[[field]]),
    logical(1)
  )))
})

test_that("E8-08 A-L1 bootstrap: unified suite replays Python-test-derived anticipation extraction", {
  skip_if_not_installed("jsonlite")

  extraction_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-a-layer1-python-tests.json"
  )

  expect_true(file.exists(extraction_path))

  extraction <- jsonlite::fromJSON(extraction_path, simplifyVector = FALSE)
  case_ids <- vapply(extraction$derived_cases, function(case) case$case_id, character(1))
  source_tests <- vapply(extraction$python_sources, function(source) source$test_name, character(1))

  expect_identical(extraction$comparison$status, "matched")
  expect_identical(extraction$comparison$exact_status, "passed")
  expect_identical(extraction$comparison$numeric_status, "not_applicable")
  expect_true(all(c("A8", "A9", "A10", "A11", "A13") %in% case_ids))
  expect_true(all(c(
    "test_sensitivity_no_anticipation",
    "test_sensitivity_no_anticipation_consistency"
  ) %in% source_tests))
  expect_equal(extraction$r_reference$api_probe$max_anticipation_tested, 3L)
})

test_that("E8-08 A-L2 bootstrap: unified suite replays common-timing Python-R sensitivity comparator", {
  skip_if_not_installed("jsonlite")

  comparator_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-a-layer2-common-timing.json"
  )

  expect_true(file.exists(comparator_path))

  comparator <- jsonlite::fromJSON(comparator_path, simplifyVector = FALSE)

  expect_identical(comparator$comparison$status, "matched")
  expect_identical(comparator$comparison$exact_status, "passed")
  expect_identical(comparator$comparison$numeric_status, "passed")
  expect_identical(comparator$python_reference$pre_period$robustness_level, "highly_sensitive")
  expect_true(isTRUE(comparator$python_reference$anticipation$anticipation_detected))
  expect_equal(comparator$python_reference$anticipation$recommended_exclusion, 1L)
})

test_that("E8-08 A-L3 bootstrap: unified suite replays smoking sensitivity oracle", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260323-qa-parity-e8-03-smoking-comparator.json"
  )

  expect_true(file.exists(oracle_path))

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  current <- compute_e8_08_smoking_summary()

  expect_true(isTRUE(oracle$comparison$exact_status == "passed"))
  expect_true(isTRUE(oracle$comparison$numeric_status == "passed"))
  expect_identical(current$robustness_level, oracle$python_reference$result$pre_period$robustness_level)
  expect_equal(
    current$sensitivity_ratio,
    oracle$python_reference$result$pre_period$sensitivity_ratio,
    tolerance = 1e-6
  )
  expect_equal(
    current$demean_att,
    oracle$python_reference$result$transformation$demean_att,
    tolerance = 1e-6
  )
})

test_that("E8-08 A-L4 bootstrap: unified suite replays Monte Carlo sensitivity oracle", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260323-qa-parity-e8-03-monte-carlo-common-timing.json"
  )
  fixture_path <- resolve_e8_08_parity_path(
    "e8_03_monte_carlo_common_timing_fixture.csv"
  )

  expect_true(file.exists(oracle_path))
  expect_true(file.exists(fixture_path))

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  fixture <- read.csv(fixture_path, stringsAsFactors = FALSE)
  current <- compute_e8_08_monte_carlo_summary(fixture)

  expect_true(isTRUE(oracle$comparison$status == "passed"))
  expect_true(isTRUE(oracle$comparison$numeric_status == "passed"))
  expect_equal(current$mean_att, oracle$python_reference$summary$mean_att, tolerance = 1e-6)
  expect_equal(current$mean_se, oracle$python_reference$summary$mean_se, tolerance = 1e-4)
  expect_equal(current$coverage, oracle$python_reference$summary$coverage, tolerance = 1e-12)
})

test_that("E8-08 A-L0 bootstrap: unified suite replays sensitivity README/example contract", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-layer0-readme-example.json"
  )

  expect_true(file.exists(oracle_path))

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  current <- compute_e8_08_readme_contract()

  expect_identical(oracle$layer, "layer_0")
  expect_identical(oracle$parity_expectation, "readme-example-contract")
  expect_identical(
    current$warnings,
    unname(unlist(oracle$example$expected_public_contract$warnings, use.names = FALSE))
  )

  expect_equal(
    current$pre_period$sensitivity_ratio,
    oracle$example$expected_public_contract$pre_period$sensitivity_ratio,
    tolerance = 1e-6
  )
  expect_identical(
    current$pre_period$robustness_level,
    oracle$example$expected_public_contract$pre_period$robustness_level
  )
  expect_identical(
    current$pre_period$is_robust,
    isTRUE(oracle$example$expected_public_contract$pre_period$is_robust)
  )

  expect_identical(
    current$anticipation$anticipation_detected,
    isTRUE(oracle$example$expected_public_contract$anticipation$anticipation_detected)
  )
  expect_equal(
    current$anticipation$recommended_exclusion,
    oracle$example$expected_public_contract$anticipation$recommended_exclusion
  )
  expect_identical(
    current$anticipation$detection_method,
    oracle$example$expected_public_contract$anticipation$detection_method
  )

  expect_true(
    is.character(current$comprehensive$overall_assessment) &&
      nchar(current$comprehensive$overall_assessment) > 0L,
    info = "overall_assessment must be a non-empty string"
  )
  expect_true(
    is.character(current$comprehensive$recommendations) &&
      all(nchar(current$comprehensive$recommendations) > 0L),
    info = "recommendations must be non-empty character"
  )
  expect_identical(
    current$comprehensive$estimator_is_null,
    isTRUE(oracle$example$expected_public_contract$comprehensive$estimator_is_null)
  )

  expected_transform <- oracle$example$expected_public_contract$comprehensive$transformation
  expect_equal(current$comprehensive$transformation$demean_att, expected_transform$demean_att, tolerance = 1e-6)
  expect_equal(current$comprehensive$transformation$demean_se, expected_transform$demean_se, tolerance = 1e-4)
  expect_equal(current$comprehensive$transformation$detrend_att, expected_transform$detrend_att, tolerance = 1e-6)
  expect_equal(current$comprehensive$transformation$detrend_se, expected_transform$detrend_se, tolerance = 1e-4)
  expect_equal(current$comprehensive$transformation$difference, expected_transform$difference, tolerance = 1e-6)
  expect_equal(current$comprehensive$transformation$rel_diff, expected_transform$rel_diff, tolerance = 1e-6)

  print_fragments <- unlist(
    oracle$example$expected_public_contract$comprehensive$print_fragments,
    use.names = FALSE
  )
  summary_fragments <- unlist(
    oracle$example$expected_public_contract$comprehensive$summary_fragments,
    use.names = FALSE
  )
  current_print <- paste(current$comprehensive$print_output, collapse = "\n")
  current_summary <- paste(current$comprehensive$summary_output, collapse = "\n")

  expect_true(all(vapply(print_fragments, grepl, logical(1), x = current_print, fixed = TRUE)))
  expect_true(all(vapply(summary_fragments, grepl, logical(1), x = current_summary, fixed = TRUE)))
})

test_that("E8-08 B1 helper: clustering stats use observation counts and reliability arithmetic", {
  oracle <- read_e8_08_required_parity_json(
    "20260327-qa-parity-e8-08-task9-numeric-oracle.json"
  )
  if (is.null(oracle)) {
    return(invisible(NULL))
  }

  hierarchy <- make_e8_08_clustering_hierarchy_panel()

  state_stats <- lwdid:::.analyze_cluster_var(
    hierarchy,
    ivar = "id",
    cluster_var = "state",
    gvar = "cohort",
    d = NULL
  )

  expected_sizes <- stats::setNames(rep(6L, 6L), sprintf("s%02d", 1:6))
  expected_reliability <- 0.5 * min(state_stats$n_clusters / 50, 1) +
    0.3 * 1 +
    0.2 * 1

  expect_equal(state_stats$cluster_sizes, expected_sizes)
  expect_equal(state_stats$mean_size, 6)
  expect_equal(state_stats$median_size, 6)
  expect_equal(state_stats$balance_score, 1, tolerance = 1e-12)
  expect_equal(state_stats$cv_score, 1, tolerance = 1e-12)
  expect_equal(
    state_stats$reliability_score,
    expected_reliability,
    tolerance = 1e-12
  )

  reliability_case <- oracle$cases$clustering_reliability_formula
  oracle_sizes <- stats::setNames(
    as.integer(unlist(reliability_case$manual_expected$cluster_sizes, use.names = FALSE)),
    names(reliability_case$manual_expected$cluster_sizes)
  )

  expect_equal(state_stats$cluster_sizes, oracle_sizes)
  expect_equal(
    state_stats$reliability_score,
    as.numeric(reliability_case$manual_expected$reliability_score),
    tolerance = 1e-12
  )
  expect_identical(reliability_case$comparison$status, "matched")
})

test_that("E8-08 B2 helper: nested clustering candidates are flagged as lower-than-unit", {
  hierarchy <- make_e8_08_clustering_hierarchy_panel()

  nested_stats <- lwdid:::.analyze_cluster_var(
    hierarchy,
    ivar = "id",
    cluster_var = "subcluster",
    gvar = "cohort",
    d = NULL
  )

  expect_true(nested_stats$is_nested_in_unit)
  expect_identical(nested_stats$level_relative_to_unit, "lower")
  expect_false(nested_stats$is_valid)
})

test_that("E8-08 B3 helper: clustering reliability scores stay within the unit interval", {
  hierarchy <- make_e8_08_clustering_hierarchy_panel()

  stats_list <- list(
    lwdid:::.analyze_cluster_var(hierarchy, "id", "region", "cohort", NULL),
    lwdid:::.analyze_cluster_var(hierarchy, "id", "state", "cohort", NULL),
    lwdid:::.analyze_cluster_var(hierarchy, "id", "subcluster", "cohort", NULL)
  )
  scores <- vapply(stats_list, function(x) x$reliability_score, numeric(1))

  expect_true(all(scores >= 0))
  expect_true(all(scores <= 1))
})

test_that("E8-08 B4 helper: treatment-variation level sorts cluster candidates by unique count", {
  hierarchy <- make_e8_08_clustering_hierarchy_panel()

  expect_identical(
    lwdid:::.detect_treatment_variation_level(
      hierarchy,
      ivar = "id",
      potential_cluster_vars = c("state", "region"),
      gvar = "cohort",
      d = NULL
    ),
    "state"
  )

  common_timing <- transform(
    hierarchy,
    d = as.integer(time >= 2L & cohort == 2005L)
  )

  expect_identical(
    lwdid:::.detect_treatment_variation_level(
      common_timing,
      ivar = "id",
      potential_cluster_vars = c("region", "state"),
      gvar = NULL,
      d = "d"
    ),
    "id"
  )
})

test_that("E8-08 B5 helper: clustering recommendation prioritizes treatment level before raw score", {
  immediate <- lwdid:::.generate_clustering_recommendation(
    cluster_structure = list(
      state = list(
        var_name = "state",
        n_clusters = 24L,
        reliability_score = 0.62,
        is_valid = TRUE,
        treatment_varies_within = FALSE,
        n_clusters_with_variation = 0L
      ),
      county = list(
        var_name = "county",
        n_clusters = 48L,
        reliability_score = 0.95,
        is_valid = TRUE,
        treatment_varies_within = FALSE,
        n_clusters_with_variation = 0L
      )
    ),
    treatment_level = "state"
  )

  expect_identical(immediate$recommended_var, "state")
  expect_identical(immediate$warnings, character(0))

  fallback <- lwdid:::.generate_clustering_recommendation(
    cluster_structure = list(
      state = list(
        var_name = "state",
        n_clusters = 15L,
        reliability_score = 0.92,
        is_valid = TRUE,
        treatment_varies_within = FALSE,
        n_clusters_with_variation = 0L
      ),
      county = list(
        var_name = "county",
        n_clusters = 25L,
        reliability_score = 0.55,
        is_valid = TRUE,
        treatment_varies_within = FALSE,
        n_clusters_with_variation = 0L
      )
    ),
    treatment_level = "state"
  )

  expect_identical(fallback$recommended_var, "county")
  expect_true(any(grepl("state", fallback$warnings, fixed = TRUE)))
})

test_that("E8-08 B6 api: few-cluster recommendation preserves wild-bootstrap contract", {
  skip_if_not_installed("jsonlite")

  oracle <- read_e8_08_parity_json("20260324-qa-parity-e8-04-layer3-few-cluster.json")
  case <- oracle$cases$few_cluster
  fixture <- read_e8_08_parity_csv(case$fixture_csv)

  recommendation <- recommend_clustering(
    fixture,
    ivar = case$ivar,
    tvar = case$tvar,
    potential_cluster_vars = unlist(case$potential_cluster_vars, use.names = FALSE),
    gvar = case$gvar,
    d = NULL,
    min_clusters = case$min_clusters,
    verbose = FALSE
  )

  expect_identical(
    recommendation$recommended_var,
    case$python_recommend_clustering$recommended_var
  )
  expect_identical(
    recommendation$n_clusters,
    case$python_recommend_clustering$n_clusters
  )
  expect_equal(
    recommendation$confidence,
    case$python_recommend_clustering$confidence,
    tolerance = 1e-12
  )
  expect_true(isTRUE(recommendation$use_wild_bootstrap))
  expect_identical(
    recommendation$wild_bootstrap_reason,
    case$python_recommend_clustering$wild_bootstrap_reason
  )
  expect_identical(
    recommendation$warnings,
    unlist(case$python_recommend_clustering$warnings, use.names = FALSE)
  )
})

test_that("E8-08 B7 api: invalid clustering candidates return an explicit no-valid-options contract", {
  hierarchy <- make_e8_08_clustering_hierarchy_panel()

  diagnosis <- diagnose_clustering(
    hierarchy,
    ivar = "id",
    potential_cluster_vars = c("subcluster"),
    gvar = "cohort",
    verbose = FALSE
  )
  recommendation <- recommend_clustering(
    hierarchy,
    ivar = "id",
    tvar = "time",
    potential_cluster_vars = c("subcluster"),
    gvar = "cohort",
    verbose = FALSE
  )

  expect_null(diagnosis$recommended_cluster_var)
  expect_identical(
    diagnosis$recommendation_reason,
    "No valid clustering options available."
  )
  expect_identical(
    diagnosis$warnings,
    "All potential cluster variables are invalid (nested within units or < 2 clusters)."
  )

  expect_null(recommendation$recommended_var)
  expect_identical(recommendation$reasons, "No valid clustering options available.")
  expect_identical(recommendation$warnings, diagnosis$warnings)
  expect_false(isTRUE(recommendation$use_wild_bootstrap))
})

test_that("E8-08 B8 api: delegated clustering recommendations cap alternatives at two ranked options", {
  hierarchy <- make_e8_08_clustering_hierarchy_panel(
    n_regions = 2L,
    states_per_region = 12L,
    units_per_state = 2L,
    periods = 2L
  )

  state_index <- as.integer(sub("^s", "", hierarchy$state))
  hierarchy$division <- sprintf("d%02d", ceiling(state_index / 6))
  hierarchy$zone <- sprintf("z%02d", ceiling(state_index / 4))
  hierarchy$block <- sprintf("b%02d", ceiling(state_index / 3))

  recommendation <- recommend_clustering(
    hierarchy,
    ivar = "id",
    tvar = "time",
    potential_cluster_vars = c("state", "block", "zone", "division"),
    gvar = "cohort",
    verbose = FALSE
  )

  expect_identical(recommendation$recommended_var, "state")
  expect_length(recommendation$alternatives, 2L)
  expect_identical(
    vapply(recommendation$alternatives, function(option) option$var, character(1)),
    c("block", "zone")
  )
})

test_that("E8-08 B-L3 bootstrap: unified suite replays hierarchical clustering parity contract", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-b-layer3-hierarchical.json"
  )

  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  current <- compute_e8_08_clustering_hierarchical_summary()

  expect_identical(oracle$comparison$status, "matched")
  expect_identical(oracle$comparison$exact_status, "passed")
  expect_identical(oracle$comparison$numeric_status, "passed")
  expect_identical(
    current$diagnosis$recommended_cluster_var,
    oracle$r_current$diagnosis$recommended_cluster_var
  )
  expect_identical(
    current$recommendation$recommended_var,
    oracle$r_current$recommendation$recommended_var
  )
  expect_identical(
    current$consistency$is_consistent,
    isTRUE(oracle$r_current$consistency$is_consistent)
  )
})

test_that("E8-08 B-L5 bootstrap: unified suite freezes current clustering release-blocker evidence", {
  skip_if_not_installed("jsonlite")

  oracle <- read_e8_08_parity_json("20260324-qa-parity-e8-04-layer5-release-regression.json")

  expect_identical(oracle$story, "story-E8-04")
  expect_identical(oracle$layer, "layer_5")
  expect_identical(oracle$exact_status, "passed")
  expect_identical(oracle$legacy_hardening_probe$status, "passed")
  expect_identical(oracle$visualization_export_suite$status, "passed")
  expect_true(isTRUE(oracle$direct_probes$by_cohort_export$literal_match))
  expect_true(isTRUE(oracle$direct_probes$ri_summary$valid_fraction_present))
  expect_identical(
    oracle$blocker_boundary,
    "targeted-release-contracts-cleared-full-check-pending"
  )
  expect_identical(oracle$remaining_gap, "none")
})

test_that("E8-08 D-L2 bootstrap: unified suite replays public trend diagnostics and recommendation parity", {
  skip_if_not_installed("jsonlite")

  oracle <- read_e8_08_parity_json("20260326-qa-parity-e8-08-section-d-layer2-public-diagnostics.json")
  wrapper_case <- oracle$cases$heterogeneous_trends_recommendation_public_api
  fixture <- read_e8_08_parity_csv(wrapper_case$fixture_csv)
  heterogeneity <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )
  recommendation <- lwdid_recommend_transformation(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    run_all_diagnostics = TRUE,
    verbose = FALSE
  )

  expected_heterogeneity <- wrapper_case$r_current_behavior$heterogeneity
  expected_recommendation <- wrapper_case$r_current_behavior$recommendation

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.5")
  expect_identical(
    oracle$source_oracle,
    "20260325-qa-parity-e8-06-trend-public-diagnostics.json"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(
    heterogeneity$has_heterogeneous_trends,
    isTRUE(expected_heterogeneity$has_heterogeneous_trends)
  )
  expect_identical(heterogeneity$recommendation, expected_heterogeneity$recommendation)
  expect_equal(
    heterogeneity$recommendation_confidence,
    as.numeric(expected_heterogeneity$recommendation_confidence),
    tolerance = 1e-10
  )
  expect_equal(
    heterogeneity$trend_heterogeneity_test$f_stat,
    as.numeric(expected_heterogeneity$trend_heterogeneity_test$f_stat),
    tolerance = 1e-10
  )
  expect_equal(
    heterogeneity$trend_heterogeneity_test$pvalue,
    as.numeric(expected_heterogeneity$trend_heterogeneity_test$pvalue),
    tolerance = 1e-10
  )
  expect_equal(
    unname(as.integer(c(
      heterogeneity$trend_heterogeneity_test$df_num,
      heterogeneity$trend_heterogeneity_test$df_den
    ))),
    as.integer(unlist(expected_heterogeneity$trend_heterogeneity_test[c("df_num", "df_den")]))
  )
  expect_equal(
    length(heterogeneity$trend_by_cohort),
    length(expected_heterogeneity$trend_by_cohort)
  )
  expect_equal(
    length(heterogeneity$trend_differences),
    length(expected_heterogeneity$trend_differences)
  )

  expect_identical(
    recommendation$recommended_method,
    expected_recommendation$recommended_method
  )
  expect_equal(
    recommendation$confidence,
    as.numeric(expected_recommendation$confidence),
    tolerance = 1e-10
  )
  expect_identical(
    recommendation$confidence_level,
    expected_recommendation$confidence_level
  )
  expect_identical(
    recommendation$alternative_method,
    expected_recommendation$alternative_method
  )
  expect_identical(
    recommendation$alternative_reason,
    expected_recommendation$alternative_reason
  )
  expect_identical(
    unname(as.character(recommendation$reasons)),
    unname(as.character(unlist(expected_recommendation$reasons, use.names = FALSE)))
  )
  expect_identical(
    unname(as.character(recommendation$warnings)),
    unname(as.character(unlist(expected_recommendation$warnings, use.names = FALSE)))
  )
  expect_equal(
    unname(as.numeric(recommendation$scores)),
    as.numeric(unlist(expected_recommendation$scores, use.names = FALSE)),
    tolerance = 1e-10
  )
  expect_identical(
    recommendation$has_seasonal_pattern,
    isTRUE(expected_recommendation$has_seasonal_pattern)
  )
  expect_identical(
    recommendation$is_balanced_panel,
    isTRUE(expected_recommendation$is_balanced_panel)
  )
})

test_that("E8-08 D11 bootstrap: unified suite freezes Welch-Satterthwaite df parity on shared trend-difference fixture", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-story-worker-e8-08-section-d-layer3-welch-df.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$heterogeneous_trends_welch_df_public_api
  fixture <- read_e8_08_parity_csv(wrapper_case$fixture_csv)
  heterogeneity <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )

  actual_differences <- lapply(heterogeneity$trend_differences, function(diff) {
    list(
      cohort_1 = as.integer(diff$cohort_1),
      cohort_2 = as.integer(diff$cohort_2),
      slope_diff = as.numeric(diff$slope_diff),
      slope_diff_se = as.numeric(diff$slope_diff_se),
      df = as.integer(diff$df)
    )
  })
  expected_differences <- wrapper_case$r_current_behavior$trend_differences
  numeric_oracle <- read_e8_08_required_parity_json(
    "20260327-qa-parity-e8-08-task9-numeric-oracle.json"
  )
  if (is.null(numeric_oracle)) {
    return(invisible(NULL))
  }

  welch_case <- numeric_oracle$cases$welch_df_formula
  all_trends <- heterogeneity$trend_by_cohort
  if (!is.null(heterogeneity$control_group_trend)) {
    all_trends <- c(all_trends, list(heterogeneity$control_group_trend))
  }

  trend_lookup <- setNames(
    all_trends,
    vapply(
      all_trends,
      function(trend) as.character(trend$cohort),
      character(1)
    )
  )
  manual_pairs <- lapply(welch_case$manual_pairs, function(pair) {
    trend_1 <- trend_lookup[[as.character(pair$cohort_1)]]
    trend_2 <- trend_lookup[[as.character(pair$cohort_2)]]
    slope_se_1 <- as.numeric(trend_1$slope_se)
    slope_se_2 <- as.numeric(trend_2$slope_se)
    numerator <- (slope_se_1^2 + slope_se_2^2)^2
    denominator <- (slope_se_1^4 / max(as.integer(trend_1$n_units) - 1L, 1L)) +
      (slope_se_2^4 / max(as.integer(trend_2$n_units) - 1L, 1L))

    list(
      cohort_1 = as.integer(pair$cohort_1),
      cohort_2 = as.integer(pair$cohort_2),
      manual_df = if (denominator > 0) {
        max(as.integer(floor(numerator / denominator)), 1L)
      } else {
        1L
      }
    )
  })

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.5")
  expect_identical(
    oracle$source_oracle,
    "20260325-story-worker-e8-06-trend-public-welch-df.json"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")
  expect_equal(length(actual_differences), length(expected_differences))

  for (i in seq_along(expected_differences)) {
    expect_identical(
      actual_differences[[i]]$cohort_1,
      as.integer(expected_differences[[i]]$cohort_1)
    )
    expect_identical(
      actual_differences[[i]]$cohort_2,
      as.integer(expected_differences[[i]]$cohort_2)
    )
    expect_equal(
      actual_differences[[i]]$slope_diff,
      as.numeric(expected_differences[[i]]$slope_diff),
      tolerance = 1e-10
    )
    expect_equal(
      actual_differences[[i]]$slope_diff_se,
      as.numeric(expected_differences[[i]]$slope_diff_se),
      tolerance = 1e-10
    )
    expect_identical(
      actual_differences[[i]]$df,
      as.integer(expected_differences[[i]]$df)
    )
  }

  expect_equal(
    vapply(manual_pairs, `[[`, integer(1), "manual_df"),
    vapply(actual_differences, `[[`, integer(1), "df")
  )
  expect_equal(
    vapply(manual_pairs, `[[`, integer(1), "manual_df"),
    vapply(welch_case$manual_pairs, function(pair) as.integer(pair$manual_df), integer(1))
  )
  expect_identical(welch_case$comparison$status, "matched")
})

test_that("E8-08 D7 bootstrap: unified suite freezes seasonal recommendation parity on a shared synthetic fixture", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-d-layer2-seasonal-recommendation.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$recommend_transformation_seasonal_public_api
  fixture <- read_e8_08_parity_csv(wrapper_case$fixture_csv)
  recommendation <- lwdid_recommend_transformation(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    run_all_diagnostics = FALSE,
    verbose = FALSE
  )

  expected <- wrapper_case$r_current_behavior
  actual_scores <- unname(as.numeric(recommendation$scores))
  expected_scores <- as.numeric(unlist(expected$scores, use.names = FALSE))

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.5")
  expect_identical(
    oracle$source_oracle,
    "20260325-qa-parity-e8-06-trend-public-seasonal.json"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(
    recommendation$recommended_method,
    expected$recommended_method
  )
  expect_equal(
    recommendation$confidence,
    as.numeric(expected$confidence),
    tolerance = 1e-10
  )
  expect_identical(
    recommendation$confidence_level,
    expected$confidence_level
  )
  expect_identical(
    recommendation$alternative_method,
    expected$alternative_method
  )
  expect_identical(
    recommendation$alternative_reason,
    expected$alternative_reason
  )
  expect_identical(
    recommendation$has_seasonal_pattern,
    isTRUE(expected$has_seasonal_pattern)
  )
  expect_identical(
    recommendation$is_balanced_panel,
    isTRUE(expected$is_balanced_panel)
  )
  expect_identical(
    unname(as.character(recommendation$reasons)),
    unname(as.character(unlist(expected$reasons, use.names = FALSE)))
  )
  expect_identical(
    unname(as.character(recommendation$warnings)),
    unname(as.character(unlist(expected$warnings, use.names = FALSE)))
  )
  expect_equal(actual_scores, expected_scores, tolerance = 1e-10)
  expect_true(all(is.finite(actual_scores)))
  expect_true(all(actual_scores >= 0 & actual_scores <= 1))
})

test_that("E8-08 D15 bootstrap: unified suite freezes seasonal-threshold helper parity and keeps cohort-trend plotting live", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("ggplot2")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-d-layer1-seasonal-helper.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  short_case <- oracle$cases$short_panel_no_seasonality
  seasonal_case <- oracle$cases$quarterly_pattern_seasonality
  short_fixture <- read_e8_08_parity_csv(short_case$fixture_csv)
  seasonal_fixture <- read_e8_08_parity_csv(seasonal_case$fixture_csv)
  plot_fixture <- read_e8_08_parity_csv(oracle$plot_fixture_csv)

  short_result <- lwdid:::.detect_seasonal_patterns(
    short_fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    threshold = as.numeric(short_case$threshold)
  )
  seasonal_result <- lwdid:::.detect_seasonal_patterns(
    seasonal_fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    threshold = as.numeric(seasonal_case$threshold)
  )
  trend_plot <- plot_cohort_trends(
    plot_fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    confidence_bands = FALSE
  )

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.5")
  expect_identical(
    oracle$source_oracle,
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )
  expect_identical(short_result, isTRUE(short_case$expected_has_seasonal))
  expect_identical(
    seasonal_result,
    isTRUE(seasonal_case$expected_has_seasonal)
  )
  expect_s3_class(trend_plot, "ggplot")
})

test_that("E8-08 E1 bootstrap: unified suite freezes return_diagnostics attachment on built-in data", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-e-return-diagnostics.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$return_diagnostics_mainline_surface
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_return_diagnostics_surface_summary()

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.3"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(
    current$castle$diagnostics_names,
    unname(as.character(expected$castle$diagnostics_names))
  )
  expect_identical(current$castle$selection_class, expected$castle$selection_class)
  expect_identical(current$castle$trends_class, expected$castle$trends_class)
  expect_identical(
    current$castle$parallel_trends_class,
    expected$castle$parallel_trends_class
  )
  expect_identical(current$castle$clustering_class, expected$castle$clustering_class)
  expect_identical(
    current$castle$controls_tier_is_null,
    isTRUE(expected$castle$controls_tier_is_null)
  )

  expect_identical(
    current$smoking$diagnostics_names,
    unname(as.character(expected$smoking$diagnostics_names))
  )
  expect_identical(
    current$smoking$controls_tier,
    as.character(expected$smoking$controls_tier)
  )
  expect_identical(current$smoking$selection_class, expected$smoking$selection_class)
  expect_identical(current$smoking$trends_class, expected$smoking$trends_class)
  expect_identical(
    current$smoking$parallel_trends_class,
    expected$smoking$parallel_trends_class
  )
  expect_identical(
    current$smoking$clustering_is_null,
    isTRUE(expected$smoking$clustering_is_null)
  )
  expect_identical(
    current$smoking$captured_messages,
    unname(as.character(expected$smoking$captured_messages))
  )
  expect_identical(
    current$smoking$captured_output,
    unname(as.character(expected$smoking$captured_output))
  )

  fatal_oracle_path <- resolve_e8_08_parity_path(
    "20260327-story-worker-e8-08-section-e-fatal-guards.json"
  )
  expect_true(file.exists(fatal_oracle_path))
  if (!file.exists(fatal_oracle_path)) {
    return(invisible(NULL))
  }

  fatal_oracle <- jsonlite::fromJSON(fatal_oracle_path, simplifyVector = FALSE)
  fatal_case <- fatal_oracle$cases$fatal_guard_surface
  fatal_expected <- fatal_case$r_current_behavior
  fatal_current <- compute_e8_08_fatal_guard_summary()

  expect_identical(fatal_oracle$story, "story-E8-08")
  expect_identical(fatal_oracle$task, "E8-08.8")
  expect_identical(fatal_case$comparison$status, "matched")
  expect_identical(fatal_case$comparison$exact_status, "passed")
  expect_identical(fatal_case$comparison$numeric_status, "passed")

  expect_identical(
    fatal_current$fatal_001$focal_cohort,
    as.integer(fatal_expected$fatal_001$focal_cohort)
  )
  expect_identical(
    fatal_current$fatal_001$pre_boundary_period,
    as.integer(fatal_expected$fatal_001$pre_boundary_period)
  )
  expect_identical(
    fatal_current$fatal_001$boundary_period,
    as.integer(fatal_expected$fatal_001$boundary_period)
  )
  expect_identical(
    fatal_current$fatal_001$pre_boundary_control_ids,
    as.integer(unlist(
      fatal_expected$fatal_001$pre_boundary_control_ids,
      use.names = FALSE
    ))
  )
  expect_identical(
    fatal_current$fatal_001$boundary_control_ids,
    as.integer(unlist(
      fatal_expected$fatal_001$boundary_control_ids,
      use.names = FALSE
    ))
  )
  expect_identical(
    fatal_current$fatal_001$pre_boundary_n_control,
    as.integer(fatal_expected$fatal_001$pre_boundary_n_control)
  )
  expect_identical(
    fatal_current$fatal_001$boundary_n_control,
    as.integer(fatal_expected$fatal_001$boundary_n_control)
  )
  expect_identical(
    fatal_current$fatal_001$pre_boundary_matches_expected,
    isTRUE(fatal_expected$fatal_001$pre_boundary_matches_expected)
  )
  expect_identical(
    fatal_current$fatal_001$boundary_matches_expected,
    isTRUE(fatal_expected$fatal_001$boundary_matches_expected)
  )

  expect_identical(
    fatal_current$fatal_004$requested_control_group,
    as.character(fatal_expected$fatal_004$requested_control_group)
  )
  expect_identical(
    fatal_current$fatal_004$used_control_group,
    as.character(fatal_expected$fatal_004$used_control_group)
  )
  expect_identical(
    fatal_current$fatal_004$n_never_treated,
    as.integer(fatal_expected$fatal_004$n_never_treated)
  )
  expect_identical(
    fatal_current$fatal_004$n_control_summary,
    as.integer(fatal_expected$fatal_004$n_control_summary)
  )
  expect_identical(
    fatal_current$fatal_004$overall_n_control,
    as.integer(fatal_expected$fatal_004$overall_n_control)
  )
  expect_identical(
    fatal_current$fatal_004$cohort_n_control,
    as.integer(unlist(
      fatal_expected$fatal_004$cohort_n_control,
      use.names = FALSE
    ))
  )
  expect_identical(
    fatal_current$fatal_004$warning_classes,
    unname(as.character(unlist(
      fatal_expected$fatal_004$warning_classes,
      use.names = FALSE
    )))
  )
  expect_identical(
    fatal_current$fatal_004$warning_messages,
    unname(as.character(unlist(
      fatal_expected$fatal_004$warning_messages,
      use.names = FALSE
    )))
  )
  expect_identical(
    fatal_current$fatal_004$switch_warning_detected,
    isTRUE(fatal_expected$fatal_004$switch_warning_detected)
  )
  expect_identical(
    fatal_current$fatal_004$nt_only_confirmed,
    isTRUE(fatal_expected$fatal_004$nt_only_confirmed)
  )
})

test_that("E8-08 E2 bootstrap: unified suite freezes mainline diagnostics invariance on built-in data", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-story-worker-e8-08-section-e-mainline-invariance.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$mainline_diagnostics_invariance
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_diagnostics_invariance_summary()

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.3"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(current$castle$diagnosed$att, current$castle$base$att)
  expect_identical(current$castle$diagnosed$se_att, current$castle$base$se_att)
  expect_identical(
    current$castle$diagnosed$df_inference,
    current$castle$base$df_inference
  )
  expect_identical(current$smoking$diagnosed$att, current$smoking$base$att)
  expect_identical(current$smoking$diagnosed$se_att, current$smoking$base$se_att)
  expect_identical(
    current$smoking$diagnosed$df_inference,
    current$smoking$base$df_inference
  )
  expect_length(current$smoking$diagnosed$captured_messages, 0L)
  expect_length(current$smoking$diagnosed$captured_output, 0L)

  expect_equal(
    current$castle$base$att,
    as.numeric(expected$castle$base$att),
    tolerance = 1e-12
  )
  expect_equal(
    current$castle$base$se_att,
    as.numeric(expected$castle$base$se_att),
    tolerance = 1e-12
  )
  expect_identical(
    current$castle$base$df_inference,
    as.integer(expected$castle$base$df_inference)
  )
  expect_equal(
    current$castle$diagnosed$att,
    as.numeric(expected$castle$diagnosed$att),
    tolerance = 1e-12
  )
  expect_equal(
    current$castle$diagnosed$se_att,
    as.numeric(expected$castle$diagnosed$se_att),
    tolerance = 1e-12
  )
  expect_identical(
    current$castle$diagnosed$df_inference,
    as.integer(expected$castle$diagnosed$df_inference)
  )
  expect_identical(
    current$castle$diagnosed$diagnostics$selection_class,
    expected$castle$diagnosed$diagnostics$selection_class
  )
  expect_identical(
    current$castle$diagnosed$diagnostics$trends_class,
    expected$castle$diagnosed$diagnostics$trends_class
  )
  expect_identical(
    current$castle$diagnosed$diagnostics$clustering_class,
    expected$castle$diagnosed$diagnostics$clustering_class
  )
  expect_identical(
    current$castle$diagnosed$diagnostics$controls_tier,
    if (is.null(expected$castle$diagnosed$diagnostics$controls_tier)) {
      NA_character_
    } else {
      as.character(expected$castle$diagnosed$diagnostics$controls_tier)
    }
  )

  expect_equal(
    current$smoking$base$att,
    as.numeric(expected$smoking$base$att),
    tolerance = 1e-12
  )
  expect_equal(
    current$smoking$base$se_att,
    as.numeric(expected$smoking$base$se_att),
    tolerance = 1e-12
  )
  expect_identical(
    current$smoking$base$df_inference,
    as.integer(expected$smoking$base$df_inference)
  )
  expect_equal(
    current$smoking$diagnosed$att,
    as.numeric(expected$smoking$diagnosed$att),
    tolerance = 1e-12
  )
  expect_equal(
    current$smoking$diagnosed$se_att,
    as.numeric(expected$smoking$diagnosed$se_att),
    tolerance = 1e-12
  )
  expect_identical(
    current$smoking$diagnosed$df_inference,
    as.integer(expected$smoking$diagnosed$df_inference)
  )
  expect_identical(
    current$smoking$diagnosed$diagnostics$selection_class,
    expected$smoking$diagnosed$diagnostics$selection_class
  )
  expect_identical(
    current$smoking$diagnosed$diagnostics$trends_class,
    expected$smoking$diagnosed$diagnostics$trends_class
  )
  expect_identical(
    current$smoking$diagnosed$diagnostics$clustering_is_null,
    isTRUE(expected$smoking$diagnosed$diagnostics$clustering_is_null)
  )
  expect_identical(
    current$smoking$diagnosed$diagnostics$controls_tier,
    if (is.null(expected$smoking$diagnosed$diagnostics$controls_tier)) {
      NA_character_
    } else {
      as.character(expected$smoking$diagnosed$diagnostics$controls_tier)
    }
  )
  expect_identical(
    current$smoking$diagnosed$captured_messages,
    unname(as.character(expected$smoking$diagnosed$captured_messages))
  )
  expect_identical(
    current$smoking$diagnosed$captured_output,
    unname(as.character(expected$smoking$diagnosed$captured_output))
  )
})

test_that("E8-08 E3 bootstrap: unified suite freezes get_diagnostics extraction on built-in data", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-story-worker-e8-08-section-e-get-diagnostics.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$get_diagnostics_mainline_surface
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_get_diagnostics_summary()

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.1"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(current$castle$all_names, unname(as.character(expected$castle$all_names)))
  expect_identical(
    current$castle$selection_class,
    expected$castle$selection_class
  )
  expect_identical(current$castle$trends_class, expected$castle$trends_class)
  expect_identical(
    current$castle$clustering_class,
    expected$castle$clustering_class
  )
  expect_identical(
    current$castle$sensitivity_is_null,
    isTRUE(expected$castle$sensitivity_is_null)
  )

  expect_identical(current$smoking$all_names, unname(as.character(expected$smoking$all_names)))
  expect_identical(
    current$smoking$selection_class,
    expected$smoking$selection_class
  )
  expect_identical(current$smoking$trends_class, expected$smoking$trends_class)
  expect_identical(
    current$smoking$clustering_is_null,
    isTRUE(expected$smoking$clustering_is_null)
  )
  expect_identical(
    current$smoking$sensitivity_is_null,
    isTRUE(expected$smoking$sensitivity_is_null)
  )
  expect_identical(
    current$smoking$captured_messages,
    unname(as.character(expected$smoking$captured_messages))
  )
  expect_identical(
    current$smoking$captured_output,
    unname(as.character(expected$smoking$captured_output))
  )
})

test_that("E8-08 E4 bootstrap: unified suite freezes the default get_diagnostics message contract when diagnostics are off", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-e-default-get-diagnostics.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$default_mainline_get_diagnostics_message
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_default_get_diagnostics_summary()

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.5.4"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(
    current$diagnostics_names,
    unname(as.character(expected$diagnostics_names))
  )
  expect_identical(
    current$controls_tier,
    as.character(expected$controls_tier)
  )
  expect_identical(
    current$returned_is_null,
    isTRUE(expected$returned_is_null)
  )
  expect_identical(
    current$captured_messages,
    unname(as.character(expected$captured_messages))
  )
})

test_that("E8-08 E5 bootstrap: unified suite freezes lwdid_diagnose mainline suite summaries on built-in data", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-story-worker-e8-08-section-e-diagnose-suite.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$diagnose_mainline_suite_surface
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_diagnose_suite_summary()

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.2"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(
    current$castle$class,
    unname(as.character(expected$castle$class))
  )
  expect_identical(
    current$castle$names,
    unname(as.character(expected$castle$names))
  )
  expect_identical(
    current$castle$selection_class,
    expected$castle$selection_class
  )
  expect_identical(
    current$castle$trends_class,
    expected$castle$trends_class
  )
  expect_identical(
    current$castle$clustering_class,
    expected$castle$clustering_class
  )
  expect_identical(
    current$castle$recommended_cluster_var,
    expected$castle$recommended_cluster_var
  )
  expect_identical(
    current$castle$treatment_variation_level,
    expected$castle$treatment_variation_level
  )
  expect_identical(
    current$castle$selection_risk,
    expected$castle$selection_risk
  )
  expect_equal(
    current$castle$missing_rate_overall,
    as.numeric(expected$castle$missing_rate_overall),
    tolerance = 1e-12
  )
  expect_identical(
    current$castle$recommended_method,
    expected$castle$recommended_method
  )
  expect_identical(
    current$castle$confidence_level,
    expected$castle$confidence_level
  )

  expect_identical(
    current$smoking$class,
    unname(as.character(expected$smoking$class))
  )
  expect_identical(
    current$smoking$names,
    unname(as.character(expected$smoking$names))
  )
  expect_identical(
    current$smoking$selection_class,
    expected$smoking$selection_class
  )
  expect_identical(
    current$smoking$trends_class,
    expected$smoking$trends_class
  )
  expect_identical(
    current$smoking$clustering_is_null,
    isTRUE(expected$smoking$clustering_is_null)
  )
  expect_identical(
    current$smoking$selection_risk,
    expected$smoking$selection_risk
  )
  expect_equal(
    current$smoking$missing_rate_overall,
    as.numeric(expected$smoking$missing_rate_overall),
    tolerance = 1e-12
  )
  expect_identical(
    current$smoking$recommended_method,
    expected$smoking$recommended_method
  )
  expect_identical(
    current$smoking$confidence_level,
    expected$smoking$confidence_level
  )
})

test_that("E8-08 E6 bootstrap: unified suite freezes standalone lwdid_diagnose module-isolation contract on built-in data", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-e-diagnose-failure-isolation.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$diagnose_failure_isolation_surface
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_diagnose_failure_isolation_summary()

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.2"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(
    current$class,
    unname(as.character(expected$class))
  )
  expect_identical(
    current$names,
    unname(as.character(expected$names))
  )
  expect_identical(
    current$selection_class,
    expected$selection_class
  )
  expect_identical(
    current$trends_class,
    expected$trends_class
  )
  expect_identical(
    current$clustering_is_null,
    isTRUE(expected$clustering_is_null)
  )
  expect_identical(
    current$selection_risk,
    expected$selection_risk
  )
  expect_identical(
    current$recommended_method,
    expected$recommended_method
  )
  expect_identical(
    current$confidence_level,
    expected$confidence_level
  )
  expect_identical(
    current$warnings_have_missing_cluster,
    isTRUE(expected$warnings_have_missing_cluster)
  )
})

test_that("E8-08 E7 bootstrap: unified suite freezes lwdid_diagnose print-summary contract on built-in data", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260327-story-worker-e8-08-section-e-diagnose-print-summary.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$diagnose_suite_print_summary_surface
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_diagnose_print_summary_summary()

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.2"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_identical(
    current$castle$print_output,
    unname(as.character(expected$castle$print_output))
  )
  expect_identical(
    current$castle$summary_output,
    unname(as.character(expected$castle$summary_output))
  )
  expect_identical(
    current$smoking$print_output,
    unname(as.character(expected$smoking$print_output))
  )
  expect_identical(
    current$smoking$summary_output,
    unname(as.character(expected$smoking$summary_output))
  )
})

test_that("E8-08 E9 bootstrap: unified suite freezes custom never-treated marker propagation on standalone diagnostics", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260327-qa-parity-e8-08-section-e-never-treated-propagation.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$diagnose_never_treated_marker_surface
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_never_treated_propagation_summary()

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.2"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  expect_equal(current, expected, tolerance = 1e-12)
  expect_identical(
    current$custom_marker$selection_risk,
    current$zero_marker$selection_risk
  )
  expect_equal(
    current$custom_marker$missing_rate_overall,
    current$zero_marker$missing_rate_overall,
    tolerance = 1e-12
  )
  expect_identical(
    current$custom_marker$recommended_method,
    current$zero_marker$recommended_method
  )
  expect_identical(
    current$custom_marker$confidence_level,
    current$zero_marker$confidence_level
  )
})

test_that("E8-08 E8 bootstrap: unified suite freezes diagnostics-plot dispatcher surface on built-in data", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  oracle_path <- resolve_e8_08_parity_path(
    "20260327-qa-parity-e8-08-section-e-diagnostics-plot.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$mainline_diagnostics_plot_surface
  expected <- wrapper_case$r_current_behavior
  current <- compute_e8_08_diagnostics_plot_summary()

  compare_plot_surface <- function(current_plot, expected_plot) {
    expect_identical(current_plot$class, unname(as.character(expected_plot$class)))
    expect_identical(
      current_plot$inherits_ggplot,
      isTRUE(expected_plot$inherits_ggplot)
    )
    expect_identical(
      current_plot$inherits_patchwork,
      isTRUE(expected_plot$inherits_patchwork)
    )
    expect_identical(current_plot$title, expected_plot$title)
    expect_identical(current_plot$x_label, expected_plot$x_label)
    expect_identical(current_plot$y_label, expected_plot$y_label)
  }

  compare_panel_surface <- function(current_plot, expected_plot) {
    expect_identical(current_plot$class, unname(as.character(expected_plot$class)))
    expect_identical(
      current_plot$inherits_ggplot,
      isTRUE(expected_plot$inherits_ggplot)
    )
    expect_identical(
      current_plot$inherits_patchwork,
      isTRUE(expected_plot$inherits_patchwork)
    )
    expect_identical(current_plot$base_title, expected_plot$base_title)
    expect_identical(
      current_plot$annotation_title,
      expected_plot$annotation_title
    )
    expect_identical(
      current_plot$patch_count,
      as.integer(expected_plot$patch_count)
    )
    expect_identical(
      current_plot$patch_titles,
      unname(as.character(expected_plot$patch_titles))
    )
  }

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.6")
  expect_identical(
    oracle$source_oracle,
    "test-diagnostics-integration.R::E8-07.3"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_identical(wrapper_case$comparison$exact_status, "passed")
  expect_identical(wrapper_case$comparison$numeric_status, "passed")

  compare_plot_surface(current$smoking_selection, expected$smoking_selection)
  compare_panel_surface(current$smoking_panel, expected$smoking_panel)
  compare_panel_surface(current$castle_panel, expected$castle_panel)
})

test_that("E8-08 C16 helper: never-treated detection honors custom markers and infinities", {
  expect_equal(
    lwdid:::.is_never_treated(
      c(NA_real_, 0, Inf, -Inf, 99, 5),
      never_treated_values = c(0, 99, Inf)
    ),
    c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  )
})

test_that("E8-08 C12 helper: missing-rate surfaces use observed Y availability on the full grid", {
  panel <- make_e8_08_selection_observed_y_panel()

  rates <- lwdid:::.compute_missing_rates(panel, ivar = "id", tvar = "time", y = "y")
  by_period <- lwdid:::.compute_missing_by_period(panel, ivar = "id", tvar = "time", y = "y")
  by_cohort <- lwdid:::.compute_missing_by_cohort(
    panel,
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    y = "y",
    never_treated_values = c(0, Inf)
  )

  expect_equal(rates$n_units, 2L)
  expect_equal(rates$n_periods, 2L)
  expect_equal(rates$n_expected, 4L)
  expect_equal(rates$n_observed, 3L)
  expect_equal(rates$overall_missing_rate, 0.25, tolerance = 1e-12)
  expect_false(rates$is_balanced)

  expect_equal(unname(by_period[["1"]]), 0.5, tolerance = 1e-12)
  expect_equal(unname(by_period[["2"]]), 0, tolerance = 1e-12)
  expect_identical(names(by_cohort), "2")
  expect_equal(unname(by_cohort[["2"]]), 0.5, tolerance = 1e-12)
})

test_that("E8-08 C9 api: unit-level observed-Y stats exclude y-NA pre-periods from usability", {
  panel <- make_e8_08_selection_observed_y_panel()

  stats <- get_unit_missing_stats(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    never_treated_values = c(0, Inf)
  )

  treated <- stats[stats$unit_id == 1, ]
  never_treated <- stats[stats$unit_id == 2, ]

  expect_equal(ncol(stats), 17L)
  expect_equal(treated$n_total_periods, 2L)
  expect_equal(treated$n_observed, 1L)
  expect_equal(treated$n_missing, 1L)
  expect_equal(treated$missing_rate, 0.5, tolerance = 1e-12)
  expect_equal(treated$first_observed, 2L)
  expect_equal(treated$last_observed, 2L)
  expect_equal(treated$observation_span, 1L)
  expect_equal(treated$cohort, 2L)
  expect_true(treated$is_treated)
  expect_equal(treated$n_pre_treatment, 0L)
  expect_equal(treated$n_post_treatment, 1L)
  expect_equal(treated$pre_treatment_missing_rate, 1, tolerance = 1e-12)
  expect_equal(treated$post_treatment_missing_rate, 0, tolerance = 1e-12)
  expect_false(treated$can_use_demean)
  expect_false(treated$can_use_detrend)
  expect_identical(treated$reason_if_excluded, "No pre-treatment observations")

  expect_false(never_treated$is_treated)
  expect_true(is.na(never_treated$cohort))
})

test_that("E8-08 C5 helper: attrition analysis uses observed Y availability for dropout timing", {
  panel <- make_e8_08_selection_attrition_balance_panel()

  attrition <- lwdid:::.compute_attrition_analysis(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = "gvar",
    d = NULL,
    never_treated_values = c(99, Inf)
  )

  expect_equal(attrition$n_complete, 1L)
  expect_equal(attrition$n_partial, 3L)
  expect_equal(attrition$overall_attrition, 0.75, tolerance = 1e-12)
  expect_equal(attrition$early_dropout_rate, 0.5, tolerance = 1e-12)
  expect_equal(attrition$late_entry_rate, 0.5, tolerance = 1e-12)
  expect_equal(attrition$dropout_before_treatment, 1L)
  expect_equal(attrition$dropout_after_treatment, 1L)
  expect_equal(attrition$treatment_related_attrition, 1, tolerance = 1e-12)
})

test_that("E8-08 C11 helper: balance statistics use observed-Y usability and min-max ratio", {
  panel <- make_e8_08_selection_attrition_balance_panel()

  balance <- lwdid:::.compute_balance_statistics(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = "gvar",
    d = NULL,
    never_treated_values = c(99, Inf)
  )

  expect_false(balance$is_balanced)
  expect_equal(balance$n_units, 4L)
  expect_equal(balance$n_periods, 4L)
  expect_equal(balance$min_obs_per_unit, 1L)
  expect_equal(balance$max_obs_per_unit, 4L)
  expect_equal(balance$mean_obs_per_unit, 2.5, tolerance = 1e-12)
  expect_equal(balance$std_obs_per_unit, stats::sd(c(2, 1, 4, 3)), tolerance = 1e-12)
  expect_equal(balance$balance_ratio, 0.25, tolerance = 1e-12)
  expect_equal(balance$n_treated_units, 3L)
  expect_equal(balance$units_below_demean_threshold, 1L)
  expect_equal(balance$units_below_detrend_threshold, 3L)
  expect_equal(balance$pct_usable_demean, 100 / 3 * 2, tolerance = 1e-12)
  expect_equal(balance$pct_usable_detrend, 0, tolerance = 1e-12)
})

test_that("E8-08 C1/C2/C3/C4 api: selection classification distinguishes balanced, MCAR, MAR, and MNAR panels", {
  cases <- make_e8_08_selection_missing_pattern_cases()

  balanced_diagnosis <- diagnose_selection_mechanism(
    cases$balanced,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = NULL,
    verbose = FALSE
  )
  mcar_pattern <- lwdid:::.classify_missing_pattern(
    cases$mcar,
    ivar = "id",
    tvar = "time",
    y = "y"
  )
  mar_pattern <- lwdid:::.classify_missing_pattern(
    cases$mar,
    ivar = "id",
    tvar = "time",
    y = "y",
    covariates = "x"
  )
  mnar_pattern <- lwdid:::.classify_missing_pattern(
    cases$mnar,
    ivar = "id",
    tvar = "time",
    y = "y"
  )

  expect_identical(balanced_diagnosis$missing_pattern, "MCAR")
  expect_identical(toupper(balanced_diagnosis$selection_risk), "LOW")
  expect_equal(balanced_diagnosis$selection_risk_score, 0)

  expect_identical(mcar_pattern$pattern, "MCAR")
  expect_identical(mar_pattern$pattern, "MAR")
  expect_identical(mnar_pattern$pattern, "MNAR")
})

test_that("E8-08 C6/C8 helper: risk scoring weights attrition, differential attrition, and balance factors", {
  risk <- lwdid:::.assess_selection_risk(
    missing_rates = list(
      overall_missing_rate = 0.35,
      is_balanced = FALSE
    ),
    missing_pattern = list(pattern = "MCAR"),
    attrition = list(
      overall_attrition = 0.35,
      dropout_before_treatment = 1L,
      dropout_after_treatment = 3L
    ),
    balance = list(
      is_balanced = FALSE,
      balance_ratio = 0.25
    )
  )

  expect_equal(risk$factors$missing_pattern, 0)
  expect_equal(risk$factors$attrition, 25)
  expect_equal(risk$factors$differential_attrition, 25)
  expect_equal(risk$factors$balance, 20)
  expect_equal(risk$score, 70)
  expect_identical(risk$risk, "high")
  expect_true(any(grepl("dropout after treatment", risk$warnings, ignore.case = TRUE)))
})

test_that("E8-08 C7/C10/C13/C14 api: diagnosis keeps usability schema, risk score, and warnings aligned", {
  observed_y_panel <- make_e8_08_selection_observed_y_panel()
  observed_y_stats <- get_unit_missing_stats(
    observed_y_panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    never_treated_values = c(0, Inf)
  )
  attrition_panel <- make_e8_08_selection_attrition_balance_panel()
  diagnosis <- diagnose_selection_mechanism(
    attrition_panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = "gvar",
    never_treated_values = c(99, Inf),
    verbose = FALSE
  )

  expect_identical(
    names(observed_y_stats),
    c(
      "unit_id",
      "cohort",
      "is_treated",
      "n_total_periods",
      "n_observed",
      "n_missing",
      "missing_rate",
      "first_observed",
      "last_observed",
      "observation_span",
      "n_pre_treatment",
      "n_post_treatment",
      "pre_treatment_missing_rate",
      "post_treatment_missing_rate",
      "can_use_demean",
      "can_use_detrend",
      "reason_if_excluded"
    )
  )
  expect_false(observed_y_stats$can_use_demean[[1L]])
  expect_false(observed_y_stats$can_use_detrend[[1L]])
  expect_identical(observed_y_stats$reason_if_excluded[[1L]], "No pre-treatment observations")

  expect_equal(
    diagnosis$balance_statistics$pct_usable_demean,
    100 / 3 * 2,
    tolerance = 1e-12
  )
  expect_equal(diagnosis$balance_statistics$pct_usable_detrend, 0, tolerance = 1e-12)
  expect_equal(diagnosis$selection_risk_score, 45)
  expect_identical(toupper(diagnosis$selection_risk), "MEDIUM")
  expect_length(diagnosis$warnings, 2L)
  expect_true(any(grepl("attrition", diagnosis$warnings, ignore.case = TRUE)))
  expect_true(any(grepl("balance", diagnosis$warnings, ignore.case = TRUE)))
})

test_that("E8-08 C15 api: never-treated filtering keeps 0, Inf, and NA cohorts out of treated summaries", {
  panel <- make_e8_08_selection_never_treated_filter_panel()

  diagnosis <- diagnose_selection_mechanism(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = "gvar",
    never_treated_values = c(0, Inf),
    verbose = FALSE
  )

  expect_equal(diagnosis$balance_statistics$n_treated_units, 1L)
  expect_identical(names(diagnosis$missing_rate_by_cohort), "3")
  expect_true(all(is.na(diagnosis$unit_stats$cohort[diagnosis$unit_stats$unit_id %in% 2:4])))
  expect_true(!any(diagnosis$unit_stats$is_treated[diagnosis$unit_stats$unit_id %in% 2:4]))
})

test_that("E8-08 D1 api: common-timing parallel-trends recommendation stays on demean when PT holds", {
  result <- lwdid_test_parallel_trends(
    make_e8_08_parallel_trends_common_panel(delta = 0),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "joint",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_parallel_trends")
  expect_false(result$reject_null)
  expect_identical(result$recommendation, "demean")
  expect_identical(result$gvar, "_gvar_dummy")
  expect_true(any(grepl("Assuming common timing", result$warnings, fixed = TRUE)))
})

test_that("E8-08 D2 api: common-timing placebo rejection upgrades the recommendation to detrend", {
  result <- lwdid_test_parallel_trends(
    make_e8_08_parallel_trends_common_panel(delta = 1),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "joint",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_true(result$reject_null)
  expect_lt(result$joint_pvalue, 0.05)
  expect_identical(result$recommendation, "detrend")
})

test_that("E8-08 D3 bootstrap: unified suite freezes common-timing joint F aggregation parity", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-story-worker-e8-08-section-d-layer2-common-timing-joint-f.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  wrapper_case <- oracle$cases$common_timing_placebo_joint_f_public_api
  fixture <- read_e8_08_parity_csv(wrapper_case$fixture_csv)
  result <- lwdid_test_parallel_trends(
    fixture,
    y = "Y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "joint",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  actual_atts <- vapply(result$pre_trend_estimates, `[[`, numeric(1), "att")
  actual_ses <- vapply(result$pre_trend_estimates, `[[`, numeric(1), "se")
  actual_dfs <- vapply(result$pre_trend_estimates, `[[`, integer(1), "df")
  manual_f_stat <- mean((actual_atts / actual_ses)^2)
  manual_df <- c(length(result$pre_trend_estimates), min(actual_dfs))

  expected <- wrapper_case$r_current_behavior
  numeric_oracle <- read_e8_08_required_parity_json(
    "20260327-qa-parity-e8-08-task9-numeric-oracle.json"
  )
  if (is.null(numeric_oracle)) {
    return(invisible(NULL))
  }

  joint_f_case <- numeric_oracle$cases$joint_f_formula

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.5")
  expect_identical(
    oracle$source_oracle,
    "20260325-qa-parity-e8-06-trend-public-common-timing.json"
  )
  expect_identical(wrapper_case$comparison$status, "matched")
  expect_equal(manual_f_stat, result$joint_f_stat, tolerance = 1e-12)
  expect_equal(result$joint_f_stat, as.numeric(expected$joint_f_stat), tolerance = 1e-12)
  expect_equal(result$joint_pvalue, as.numeric(expected$joint_pvalue), tolerance = 1e-12)
  expect_equal(unname(as.integer(result$joint_df)), as.integer(unlist(expected$joint_df)))
  expect_equal(unname(as.integer(manual_df)), as.integer(unlist(expected$joint_df)))
  expect_equal(
    manual_f_stat,
    as.numeric(joint_f_case$manual_expected$joint_f_stat),
    tolerance = 1e-12
  )
  expect_equal(
    as.integer(manual_df),
    as.integer(unlist(joint_f_case$manual_expected$joint_df, use.names = FALSE))
  )
  expect_identical(joint_f_case$comparison$status, "matched")
  expect_identical(result$recommendation, expected$recommendation)
  expect_identical(result$reject_null, isTRUE(expected$reject_null))
})

test_that("E8-08 D4 api: heterogeneous-trend public diagnostics keep the detrend recommendation on the shared synthetic fixture", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-d-layer2-heterogeneity-public-diagnostics.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  case <- oracle$cases$heterogeneous_trends_public_api
  fixture <- read_e8_08_parity_csv(case$fixture_csv)
  result <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.5")
  expect_identical(
    oracle$source_oracle,
    "20260325-qa-parity-e8-06-trend-public-diagnostics.json"
  )
  expect_identical(case$comparison$status, "matched")
  expect_true(result$has_heterogeneous_trends)
  expect_identical(result$recommendation, expected$recommendation)
  expect_equal(
    unname(as.numeric(result$recommendation_confidence)),
    as.numeric(expected$recommendation_confidence),
    tolerance = 1e-10
  )
  expect_equal(
    list(
      f_stat = unname(as.numeric(result$trend_heterogeneity_test$f_stat)),
      pvalue = unname(as.numeric(result$trend_heterogeneity_test$pvalue)),
      df_num = as.integer(result$trend_heterogeneity_test$df_num),
      df_den = as.integer(result$trend_heterogeneity_test$df_den),
      reject_null = isTRUE(result$trend_heterogeneity_test$reject_null)
    ),
    expected$trend_heterogeneity_test,
    tolerance = 1e-10
  )
})

test_that("E8-08 D5 api: pairwise trend differences preserve slope-difference, Welch df, and public fields", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-d-layer2-heterogeneity-public-diagnostics.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  case <- oracle$cases$heterogeneous_trends_public_api
  fixture <- read_e8_08_parity_csv(case$fixture_csv)
  result <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior$trend_differences
  actual <- lapply(result$trend_differences, function(difference) {
    list(
      cohort_1 = as.integer(difference$cohort_1),
      cohort_2 = as.integer(difference$cohort_2),
      slope_diff = as.numeric(difference$slope_diff),
      slope_diff_se = as.numeric(difference$slope_diff_se),
      pvalue = as.numeric(difference$pvalue),
      df = as.integer(difference$df)
    )
  })

  expect_true(all(c(
    "slope_1",
    "slope_2",
    "slope_diff",
    "slope_diff_se",
    "t_stat",
    "pvalue",
    "df"
  ) %in% names(result$trend_differences[[1L]])))
  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("E8-08 D6 api: transformation recommendation downgrades to demean when detrending is infeasible", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-d-layer2-detrend-infeasible.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  case <- oracle$cases$common_timing_detrend_infeasible_public_api
  fixture <- read_e8_08_parity_csv(case$fixture_csv)
  result <- lwdid_recommend_transformation(
    fixture,
    y = case$y,
    ivar = case$ivar,
    tvar = case$tvar,
    gvar = NULL,
    run_all_diagnostics = TRUE,
    verbose = FALSE
  )
  expected <- case$r_current_behavior

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.5")
  expect_identical(
    oracle$source_oracle,
    "20260326-theory-parity-e8-08-section-d-state-convergence.md"
  )
  expect_identical(case$comparison$status, "matched")
  expect_identical(case$comparison$exact_status, "passed")
  expect_identical(case$comparison$numeric_status, "passed")
  expect_s3_class(result, "lwdid_transformation_recommendation")
  expect_identical(result$recommended_method, expected$recommended_method)
  expect_equal(result$confidence, as.numeric(expected$confidence), tolerance = 1e-12)
  expect_identical(result$confidence_level, expected$confidence_level)
  expect_identical(
    is.null(result$parallel_trends_test),
    isTRUE(expected$parallel_trends_test_is_null)
  )
  expect_identical(
    is.null(result$heterogeneous_trends_diag),
    isTRUE(expected$heterogeneous_trends_diag_is_null)
  )
  expect_identical(result$n_pre_periods_min, as.integer(expected$n_pre_periods_min))
  expect_identical(result$n_pre_periods_max, as.integer(expected$n_pre_periods_max))
  expect_identical(
    isTRUE(result$has_seasonal_pattern),
    isTRUE(expected$has_seasonal_pattern)
  )
  expect_identical(
    isTRUE(result$is_balanced_panel),
    isTRUE(expected$is_balanced_panel)
  )
  expect_identical(result$alternative_method, expected$alternative_method)
  expect_identical(result$alternative_reason, expected$alternative_reason)
  expect_identical(result$warnings, expected$warnings)
  expect_identical(result$reasons, expected$reasons)
  expect_equal(
    unname(as.numeric(result$scores[c("demean", "detrend")])),
    as.numeric(unlist(expected$scores[c("demean", "detrend")], use.names = FALSE)),
    tolerance = 1e-12
  )
})

test_that("E8-08 D8 bootstrap: transformation recommendation wrapper keeps normalized scores on the unified surface", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-d-layer2-heterogeneity-public-diagnostics.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  case <- oracle$cases$recommend_transformation_public_api
  fixture <- read_e8_08_parity_csv(case$fixture_csv)
  result <- lwdid_recommend_transformation(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    run_all_diagnostics = TRUE,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  expect_identical(case$comparison$status, "matched")
  expect_identical(result$recommended_method, expected$recommended_method)
  expect_equal(unname(as.numeric(result$confidence)), as.numeric(expected$confidence), tolerance = 1e-10)
  expect_identical(result$confidence_level, expected$confidence_level)
  expect_identical(result$alternative_method, expected$alternative_method)
  expect_identical(result$alternative_reason, expected$alternative_reason)
  expect_true(all(is.finite(unname(as.numeric(result$scores)))))
  expect_true(all(unname(as.numeric(result$scores)) >= 0))
  expect_true(all(unname(as.numeric(result$scores)) <= 1))
  expect_equal(
    unname(as.numeric(result$scores[c("demean", "detrend")])),
    unname(as.numeric(unlist(expected$scores[c("demean", "detrend")]))),
    tolerance = 1e-10
  )
})

test_that("E8-08 D9 api: cohort-specific pre-trend estimates stay on the pooled-OLS public surface", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-d-layer2-heterogeneity-public-diagnostics.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  case <- oracle$cases$heterogeneous_trends_public_api
  fixture <- read_e8_08_parity_csv(case$fixture_csv)
  result <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior$trend_by_cohort
  actual <- lapply(result$trend_by_cohort, function(trend) {
    list(
      cohort = as.integer(trend$cohort),
      slope = as.numeric(trend$slope),
      slope_se = as.numeric(trend$slope_se),
      slope_pvalue = as.numeric(trend$slope_pvalue),
      n_units = as.integer(trend$n_units),
      n_pre_periods = as.integer(trend$n_pre_periods)
    )
  })

  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("E8-08 D10 api: control-group pre-trend estimates stay on the pooled-OLS public surface", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-qa-parity-e8-08-section-d-layer2-heterogeneity-public-diagnostics.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  case <- oracle$cases$heterogeneous_trends_public_api
  fixture <- read_e8_08_parity_csv(case$fixture_csv)
  result <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior$control_group_trend
  actual <- list(
    cohort = as.integer(result$control_group_trend$cohort),
    slope = as.numeric(result$control_group_trend$slope),
    slope_se = as.numeric(result$control_group_trend$slope_se),
    slope_pvalue = as.numeric(result$control_group_trend$slope_pvalue),
    n_units = as.integer(result$control_group_trend$n_units),
    n_pre_periods = as.integer(result$control_group_trend$n_pre_periods)
  )

  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("E8-08 D14 api: trend diagnostics keep cohort tags and figure placeholders on the unified surface", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_e8_08_parity_path(
    "20260326-story-worker-e8-08-section-d-layer2-common-timing-joint-f.json"
  )
  expect_true(file.exists(oracle_path))
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  common_case <- oracle$cases$common_timing_placebo_joint_f_public_api
  heterogeneity_case <- oracle$cases$heterogeneous_trends_structure_public_api
  common_fixture <- read_e8_08_parity_csv(common_case$fixture_csv)
  heterogeneity_fixture <- read_e8_08_parity_csv(heterogeneity_case$fixture_csv)

  parallel_result <- lwdid_test_parallel_trends(
    common_fixture,
    y = "Y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "joint",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )
  heterogeneity_result <- lwdid_diagnose_heterogeneous_trends(
    heterogeneity_fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )

  expected_common <- common_case$structure_contract
  expected_heterogeneity <- heterogeneity_case$structure_contract
  actual_pre_trend_fields <- names(parallel_result$pre_trend_estimates[[1L]])
  actual_pre_trend_cohorts <- vapply(
    parallel_result$pre_trend_estimates,
    `[[`,
    integer(1),
    "cohort"
  )
  actual_trend_cohorts <- vapply(
    heterogeneity_result$trend_by_cohort,
    `[[`,
    integer(1),
    "cohort"
  )

  expect_equal(
    unname(actual_pre_trend_cohorts),
    as.integer(unlist(expected_common$pre_trend_cohorts))
  )
  expect_true(all(expected_common$required_pre_trend_fields %in% actual_pre_trend_fields))
  expect_identical(
    is.null(parallel_result$figure),
    isTRUE(expected_common$parallel_trends_figure_is_null)
  )
  expect_equal(
    unname(actual_trend_cohorts),
    as.integer(unlist(expected_heterogeneity$trend_by_cohort_cohorts))
  )
  expect_identical(
    heterogeneity_result$control_group_trend$cohort,
    as.integer(expected_heterogeneity$control_group_cohort)
  )
  expect_identical(
    is.null(heterogeneity_result$figure),
    isTRUE(expected_heterogeneity$heterogeneous_trends_figure_is_null)
  )
})

test_that("E8-08 D13 api: trend diagnostics reject non-numeric never-treated markers before cohort processing", {
  panel <- make_heterogeneous_trend_panel()

  expect_error(
    lwdid_test_parallel_trends(
      panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "cohort",
      never_treated_values = "zero",
      method = "joint",
      verbose = FALSE
    ),
    class = "lwdid_invalid_parameter"
  )

  expect_error(
    lwdid_diagnose_heterogeneous_trends(
      panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "cohort",
      never_treated_values = "zero",
      verbose = FALSE
    ),
    class = "lwdid_invalid_parameter"
  )

  expect_error(
    lwdid_recommend_transformation(
      panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "cohort",
      never_treated_values = "zero",
      run_all_diagnostics = FALSE,
      verbose = FALSE
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("E8-08 D12 api: common-timing placebo estimates exclude the e=-1 anchor from the unified surface", {
  result <- lwdid_test_parallel_trends(
    make_e8_08_parallel_trends_common_panel(delta = 0),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "placebo",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_equal(
    vapply(result$pre_trend_estimates, `[[`, integer(1), "event_time"),
    c(-2L, -3L)
  )
  expect_false(any(vapply(
    result$pre_trend_estimates,
    function(estimate) identical(estimate$event_time, -1L),
    logical(1)
  )))
})

test_that("E8-08 F1 api: selection diagnostics reject empty panel inputs before observed-Y bookkeeping", {
  empty_panel <- data.frame(
    id = integer(0),
    time = integer(0),
    gvar = integer(0),
    y = numeric(0)
  )

  expect_error(
    lwdid:::.compute_missing_rates(
      empty_panel,
      ivar = "id",
      tvar = "time",
      y = "y"
    ),
    "Input data is empty",
    class = "lwdid_invalid_parameter"
  )

  expect_error(
    diagnose_selection_mechanism(
      empty_panel,
      ivar = "id",
      tvar = "time",
      y = "y",
      gvar = "gvar",
      verbose = FALSE
    ),
    "Input data is empty",
    class = "lwdid_invalid_parameter"
  )
})

test_that("E8-08 F2 api: single treated cohort keeps heterogeneity diagnostics non-rejective and omits pairwise differences", {
  result <- lwdid_diagnose_heterogeneous_trends(
    make_single_cohort_panel(),
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "cohort",
    include_control_group = TRUE,
    alpha = 0.05,
    verbose = FALSE
  )

  expect_false(result$has_heterogeneous_trends)
  expect_identical(result$recommendation, "demean")
  expect_identical(result$trend_heterogeneity_test$reject_null, FALSE)
  expect_length(result$trend_by_cohort, 1L)
  expect_identical(result$trend_by_cohort[[1L]]$cohort, 5L)
  expect_identical(result$control_group_trend$cohort, 0L)
  expect_length(result$trend_differences, 0L)
})

test_that("E8-08 F3 api: all-NA outcomes surface boundary-safe selection and trend diagnostics", {
  panel <- make_e8_08_all_na_outcome_panel()

  selection <- diagnose_selection_mechanism(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "cohort",
    verbose = FALSE
  )

  trend_warnings <- character(0)
  recommendation <- withCallingHandlers(
    lwdid_recommend_transformation(
      panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "cohort",
      run_all_diagnostics = TRUE,
      verbose = FALSE
    ),
    warning = function(w) {
      trend_warnings <<- c(trend_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(selection$missing_pattern, "UNKNOWN")
  expect_equal(selection$missing_rate_overall, 1, tolerance = 1e-12)
  expect_identical(selection$selection_risk, "high")
  expect_match(
    selection$missing_pattern_description,
    "All outcomes are missing"
  )
  expect_identical(recommendation$recommended_method, "demean")
  expect_identical(recommendation$confidence_level, "Medium")
  expect_length(recommendation$warnings, 0L)
  expect_length(trend_warnings, 0L)
  expect_true(any(grepl(
    "No valid pre-treatment estimates computed",
    recommendation$parallel_trends_test$warnings,
    fixed = TRUE
  )))
})

test_that("E8-08 F4 api: Inf and NaN outcomes do not crash selection or trend diagnostics", {
  panel <- make_e8_08_inf_nan_outcome_panel()

  selection <- diagnose_selection_mechanism(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "cohort",
    verbose = FALSE
  )
  trend_warnings <- character(0)
  recommendation <- withCallingHandlers(
    lwdid_recommend_transformation(
      panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "cohort",
      run_all_diagnostics = TRUE,
      verbose = FALSE
    ),
    warning = function(w) {
      trend_warnings <<- c(trend_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(selection$missing_pattern, "MCAR")
  expect_identical(selection$selection_risk, "medium")
  expect_equal(selection$missing_rate_overall, 2 / 24, tolerance = 1e-12)
  expect_true(any(grepl("High attrition rate", selection$warnings)))
  expect_true(any(grepl("NA/NaN/Inf in 'y'", trend_warnings, fixed = TRUE)))
  expect_identical(recommendation$recommended_method, "demean")
  expect_identical(recommendation$confidence_level, "Medium")
  expect_true(all(is.finite(unname(recommendation$scores))))
  expect_equal(sum(unname(recommendation$scores)), 1, tolerance = 1e-12)
})

test_that("E8-08 F5 api: single-period panels keep selection stable and mark detrending infeasible", {
  panel <- make_e8_08_single_period_panel()

  selection <- diagnose_selection_mechanism(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "cohort",
    verbose = FALSE
  )
  recommendation <- lwdid_recommend_transformation(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "cohort",
    run_all_diagnostics = TRUE,
    verbose = FALSE
  )
  parallel_trends <- lwdid_test_parallel_trends(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "cohort",
    verbose = FALSE
  )

  expect_identical(selection$missing_pattern, "MCAR")
  expect_equal(selection$missing_rate_overall, 0, tolerance = 1e-12)
  expect_identical(selection$selection_risk, "low")
  expect_identical(recommendation$recommended_method, "demean")
  expect_identical(recommendation$confidence_level, "High")
  expect_true(any(grepl("Detrending requires >=2 pre-treatment periods", recommendation$warnings)))
  expect_identical(parallel_trends$reject_null, FALSE)
  expect_true(any(grepl("No valid pre-treatment estimates", parallel_trends$warnings)))
})
