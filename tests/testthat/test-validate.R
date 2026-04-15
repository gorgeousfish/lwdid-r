# ============================================================================
# test-validate.R — Comprehensive tests for validate.R
# ============================================================================

# Helper: create valid Common Timing panel data
make_ct_data <- function(n_units = 10, n_periods = 8,
                         n_treated = 5, K = 4) {
  ids <- seq_len(n_units)
  years <- 2000 + seq_len(n_periods)
  df <- expand.grid(id = ids, year = years)
  df <- df[order(df$id, df$year), ]
  # Treatment: first n_treated units are treated
  df$d <- ifelse(df$id <= n_treated, 1, 0)
  # Post: periods after K
  df$post <- ifelse(df$year > (2000 + K), 1, 0)
  # Outcome: simple DGP
  set.seed(42)
  df$y <- rnorm(nrow(df)) + df$d * df$post * 2
  rownames(df) <- NULL
  df
}

# Helper: create valid Staggered panel data
make_stag_data <- function(n_units = 12, n_periods = 10,
                           cohorts = c(2005, 2007),
                           n_nt = 4) {
  ids <- seq_len(n_units)
  years <- 2000 + seq_len(n_periods)
  df <- expand.grid(id = ids, year = years)
  df <- df[order(df$id, df$year), ]
  # Assign gvar: first n_nt are NT, rest split among cohorts
  treated_ids <- (n_nt + 1):n_units
  n_per_cohort <- length(treated_ids) %/% length(cohorts)
  gvar_map <- numeric(n_units)
  gvar_map[1:n_nt] <- 0  # never-treated
  for (i in seq_along(cohorts)) {
    start_idx <- n_nt + (i - 1) * n_per_cohort + 1
    end_idx <- if (i == length(cohorts)) n_units
               else n_nt + i * n_per_cohort
    gvar_map[start_idx:end_idx] <- cohorts[i]
  }
  df$gvar <- gvar_map[df$id]
  set.seed(42)
  df$y <- rnorm(nrow(df))
  rownames(df) <- NULL
  df
}


# ============================================================================
# Group 1: Layer 1 — Reserved columns + empty data
# ============================================================================

test_that("Layer 1: empty data.frame raises lwdid_invalid_parameter", {
  df <- data.frame()
  expect_error(
    .validate_data_not_empty(df),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 1: non-data.frame raises lwdid_invalid_parameter", {
  m <- matrix(1:4, nrow = 2)
  expect_error(
    .validate_data_not_empty(m),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 1: valid data.frame passes", {
  df <- data.frame(x = 1:3)
  expect_silent(.validate_data_not_empty(df))
})

test_that("Layer 1: reserved column 'd_' raises error", {
  df <- data.frame(d_ = 1:3, y = 1:3)
  expect_error(
    .validate_reserved_columns(df),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 1: reserved column 'tindex' raises error", {
  df <- data.frame(tindex = 1:3, y = 1:3)
  expect_error(
    .validate_reserved_columns(df),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 1: multiple reserved columns in error message", {
  df <- data.frame(d_ = 1, post_ = 1, tindex = 1)
  expect_error(
    .validate_reserved_columns(df),
    regexp = "d_.*post_.*tindex|tindex.*post_.*d_",
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 1: normal data passes reserved check", {
  df <- data.frame(id = 1:3, year = 2001:2003, y = rnorm(3))
  expect_silent(.validate_reserved_columns(df))
})


# ============================================================================
# Group 2: Layer 2 — Parameter type/range validation
# ============================================================================

test_that("Layer 2: y must be string", {
  df <- data.frame(id = 1, year = 2001, y = 1)
  expect_error(
    .validate_string_param(123, "y", df),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 2: y not in data raises lwdid_missing_column", {
  df <- data.frame(id = 1, year = 2001, y = 1)
  expect_error(
    .validate_string_param("outcome", "y", df),
    class = "lwdid_missing_column"
  )
})

test_that("Layer 2: ivar not in data raises lwdid_missing_column", {
  df <- data.frame(id = 1, year = 2001, y = 1)
  expect_error(
    .validate_string_param("unit_id", "ivar", df),
    class = "lwdid_missing_column"
  )
})

test_that("Layer 2: tvar length 3 raises error", {
  df <- data.frame(a = 1, b = 2, c = 3)
  expect_error(
    .validate_tvar(c("a", "b", "c"), df),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 2: tvar length 2 with invalid quarter values", {
  df <- data.frame(year = 2001, quarter = 5)
  expect_error(
    .validate_tvar(c("year", "quarter"), df),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 2: rolling case-insensitive", {
  expect_equal(
    .validate_rolling_parameter("Demean", "year", NULL),
    "demean"
  )
  expect_equal(
    .validate_rolling_parameter("DETREND", "year", NULL),
    "detrend"
  )
  expect_equal(
    .validate_rolling_parameter("DeMeanQ", c("year", "q"), NULL),
    "demeanq"
  )
})

test_that("Layer 2: invalid rolling raises lwdid_invalid_rolling", {
  expect_error(
    .validate_rolling_parameter("invalid_method", "year", NULL),
    class = "lwdid_invalid_rolling"
  )
})


test_that("Layer 2: invalid estimator raises error", {
  expect_error(
    .validate_choice("ols", VALID_ESTIMATORS, "estimator"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 2: alpha boundaries", {
  expect_error(
    .validate_numeric_range(0, 0, 1, "alpha", exclusive = TRUE),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    .validate_numeric_range(1, 0, 1, "alpha", exclusive = TRUE),
    class = "lwdid_invalid_parameter"
  )
  expect_silent(
    .validate_numeric_range(0.5, 0, 1, "alpha", exclusive = TRUE)
  )
})

test_that("Layer 2: trim_threshold boundaries", {
  expect_error(
    .validate_numeric_range(0, 0, 0.5, "trim_threshold",
                            exclusive = TRUE),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    .validate_numeric_range(0.5, 0, 0.5, "trim_threshold",
                            exclusive = TRUE),
    class = "lwdid_invalid_parameter"
  )
  expect_silent(
    .validate_numeric_range(0.25, 0, 0.5, "trim_threshold",
                            exclusive = TRUE)
  )
})

test_that("Layer 2: n_neighbors must be positive integer", {
  expect_error(
    .validate_positive_integer(0, "n_neighbors"),
    class = "lwdid_invalid_parameter"
  )
  expect_silent(.validate_positive_integer(1, "n_neighbors"))
})

test_that("Layer 2: Q must be >= 2", {
  expect_error(
    .validate_positive_integer(1, "Q", min_val = 2L),
    class = "lwdid_invalid_parameter"
  )
  expect_silent(.validate_positive_integer(2, "Q", min_val = 2L))
  expect_silent(.validate_positive_integer(4, "Q", min_val = 2L))
})

test_that("Layer 2: vce validation", {
  expect_silent(.validate_vce(NULL))
  expect_silent(.validate_vce("robust"))
  expect_error(
    .validate_vce("invalid"),
    class = "lwdid_invalid_vce"
  )
})

test_that("Layer 2: logical validation", {
  expect_silent(.validate_logical(TRUE, "ri"))
  expect_error(
    .validate_logical("yes", "ri"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 2: balanced_panel invalid value", {
  expect_error(
    .validate_choice("invalid", VALID_BALANCED_PANEL,
                     "balanced_panel"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 2: verbose invalid value", {
  expect_error(
    .validate_choice("invalid", VALID_VERBOSE, "verbose"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 2: season_var values in range pass", {
  df <- data.frame(q = c(1, 2, 3, 4))
  expect_silent(.validate_season_var_values(df, "q", 4))
})

test_that("Layer 2: season_var out of range raises error", {
  df <- data.frame(q = c(1, 2, 5))
  expect_error(
    .validate_season_var_values(df, "q", 4),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 2: season_var NULL skips validation", {
  df <- data.frame(q = c(1, 2, 5))
  expect_silent(.validate_season_var_values(df, NULL, 4))
})


# ============================================================================
# Group 3: Layer 3 — Mode identification + binarization
# ============================================================================

test_that("Layer 3a: gvar non-NULL → staggered", {
  expect_equal(.identify_mode(NULL, NULL, "gvar_col"), "staggered")
})

test_that("Layer 3a: d+post non-NULL → common_timing", {
  expect_equal(.identify_mode("d", "post", NULL), "common_timing")
})

test_that("Layer 3a: neither → error", {
  expect_error(
    .identify_mode(NULL, NULL, NULL),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 3a: gvar+d → staggered + warning", {
  expect_warning(
    result <- .identify_mode("d", NULL, "gvar_col"),
    class = "lwdid_data"
  )
  expect_equal(result, "staggered")
})

test_that("Layer 3b: binarize_with_na preserves NA", {
  result <- .binarize_with_na(c(0, 5, NA, -1, 0.5, NaN))
  expect_equal(result, c(0L, 1L, NA_integer_, 1L, 1L, NA_integer_))
})

test_that("Layer 3b: binarize returns integer type", {
  result <- .binarize_with_na(c(0, 1, 2))
  expect_type(result, "integer")
})

test_that("Layer 3b: binarize d=0 → 0, d=5 → 1", {
  expect_equal(.binarize_with_na(0), 0L)
  expect_equal(.binarize_with_na(5), 1L)
})

test_that("Layer 3b: binarize d=-1 → 1, d=0.5 → 1", {
  expect_equal(.binarize_with_na(-1), 1L)
  expect_equal(.binarize_with_na(0.5), 1L)
})


# ============================================================================
# Group 4: Layer 4-5 — Data type + time-invariance
# ============================================================================

test_that("Layer 4: numeric y passes", {
  df <- data.frame(y = c(1.0, 2.0, 3.0))
  expect_silent(.validate_outcome_dtype(df, "y"))
})

test_that("Layer 4: logical y warns lwdid_data", {
  df <- data.frame(y = c(TRUE, FALSE, TRUE))
  expect_warning(
    .validate_outcome_dtype(df, "y"),
    class = "lwdid_data"
  )
})

test_that("Layer 4: character y raises error", {
  df <- data.frame(y = c("a", "b", "c"),
                   stringsAsFactors = FALSE)
  expect_error(
    .validate_outcome_dtype(df, "y"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 4: factor y raises error", {
  df <- data.frame(y = factor(c("a", "b")))
  expect_error(
    .validate_outcome_dtype(df, "y"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 4: controls with character raises error", {
  df <- data.frame(x1 = 1:3,
                   x2 = c("a", "b", "c"),
                   stringsAsFactors = FALSE)
  expect_error(
    .validate_controls_dtype(df, c("x1", "x2")),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 4: controls all numeric passes", {
  df <- data.frame(x1 = 1:3, x2 = 4:6)
  expect_silent(.validate_controls_dtype(df, c("x1", "x2")))
})

test_that("Layer 4: controls NULL skips", {
  df <- data.frame(x1 = 1:3)
  expect_silent(.validate_controls_dtype(df, NULL))
})


test_that("Layer 5: d time-invariant passes", {
  df <- data.frame(
    id = rep(1:3, each = 4),
    d = rep(c(1, 1, 0), each = 4)
  )
  expect_silent(
    .validate_treatment_time_invariance(df, "d", "id")
  )
})

test_that("Layer 5: d time-varying raises error", {
  df <- data.frame(
    id = rep(1:2, each = 4),
    d = c(1, 1, 0, 0, 0, 0, 0, 0)  # unit 1 varies
  )
  expect_error(
    .validate_treatment_time_invariance(df, "d", "id"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 5: floating-point precision boundary", {
  # std exactly at threshold (1e-10) should pass (> not >=)
  # Create data where within-unit std is exactly 1e-10
  # sd of c(x, x+delta) = delta/sqrt(2)
  # We want sd = 1e-10, so delta = 1e-10 * sqrt(2)
  delta_at <- 1e-10 * sqrt(2)
  df_at <- data.frame(
    id = rep(1, 2),
    d = c(0.5, 0.5 + delta_at)
  )
  # std = 1e-10, which is NOT > 1e-10, so should pass
  expect_silent(
    .validate_treatment_time_invariance(df_at, "d", "id")
  )

  # std slightly above threshold should fail
  delta_above <- 1.1e-10 * sqrt(2)
  df_above <- data.frame(
    id = rep(1, 2),
    d = c(0.5, 0.5 + delta_above)
  )
  expect_error(
    .validate_treatment_time_invariance(df_above, "d", "id"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 5: controls time-varying raises warning", {
  df <- data.frame(
    id = rep(1:2, each = 3),
    x1 = c(1, 2, 3, 4, 4, 4)  # unit 1 varies
  )
  expect_warning(
    .validate_time_invariant_controls(df, "id", "x1"),
    class = "lwdid_data"
  )
})

test_that("Layer 5: controls NULL skips", {
  df <- data.frame(id = 1:3)
  expect_silent(
    .validate_time_invariant_controls(df, "id", NULL)
  )
})


# ============================================================================
# Group 5: Layer 6 — Cross-parameter consistency
# ============================================================================

test_that("Layer 6: vce=cluster + cluster_var=NULL → error", {
  expect_error(
    .validate_cross_param_consistency(
      vce = "cluster", cluster_var = NULL, rolling = "demean",
      season_var = NULL, tvar = "year", estimator = "ra",
      controls = NULL, ps_controls = NULL, aggregate = "none",
      control_group = "not_yet_treated", graph = FALSE,
      ri = FALSE, rireps = 1000L, ri_method = "bootstrap",
      exclude_pre_periods = 0L, Q = 4L, mode = "common_timing"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 6: vce=bootstrap + cluster_var=NULL → error", {
  expect_error(
    .validate_cross_param_consistency(
      vce = "bootstrap", cluster_var = NULL, rolling = "demean",
      season_var = NULL, tvar = "year", estimator = "ra",
      controls = NULL, ps_controls = NULL, aggregate = "none",
      control_group = "not_yet_treated", graph = FALSE,
      ri = FALSE, rireps = 1000L, ri_method = "bootstrap",
      exclude_pre_periods = 0L, Q = 4L, mode = "common_timing"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 6: estimator=ra + ps_controls → warning", {
  expect_warning(
    .validate_cross_param_consistency(
      vce = NULL, cluster_var = NULL, rolling = "demean",
      season_var = NULL, tvar = "year", estimator = "ra",
      controls = NULL, ps_controls = c("x1"),
      aggregate = "none", control_group = "not_yet_treated",
      graph = FALSE, ri = FALSE, rireps = 1000L,
      ri_method = "bootstrap", exclude_pre_periods = 0L,
      Q = 4L, mode = "common_timing"
    ),
    class = "lwdid_data"
  )
})

test_that("Layer 6: estimator=psm + no controls → error", {
  expect_error(
    .validate_cross_param_consistency(
      vce = NULL, cluster_var = NULL, rolling = "demean",
      season_var = NULL, tvar = "year", estimator = "psm",
      controls = NULL, ps_controls = NULL, aggregate = "none",
      control_group = "not_yet_treated", graph = FALSE,
      ri = FALSE, rireps = 1000L, ri_method = "bootstrap",
      exclude_pre_periods = 0L, Q = 4L, mode = "common_timing"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 6: estimator=ipwra + no controls → error", {
  expect_error(
    .validate_cross_param_consistency(
      vce = NULL, cluster_var = NULL, rolling = "demean",
      season_var = NULL, tvar = "year", estimator = "ipwra",
      controls = NULL, ps_controls = c("x1"),
      aggregate = "none", control_group = "not_yet_treated",
      graph = FALSE, ri = FALSE, rireps = 1000L,
      ri_method = "bootstrap", exclude_pre_periods = 0L,
      Q = 4L, mode = "common_timing"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 6: estimator=ipw + no controls/ps_controls → error", {
  expect_error(
    .validate_cross_param_consistency(
      vce = NULL, cluster_var = NULL, rolling = "demean",
      season_var = NULL, tvar = "year", estimator = "ipw",
      controls = NULL, ps_controls = NULL, aggregate = "none",
      control_group = "not_yet_treated", graph = FALSE,
      ri = FALSE, rireps = 1000L, ri_method = "bootstrap",
      exclude_pre_periods = 0L, Q = 4L, mode = "common_timing"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 6: vce=robust + cluster_var → warning", {
  expect_warning(
    .validate_cross_param_consistency(
      vce = "robust", cluster_var = "cl", rolling = "demean",
      season_var = NULL, tvar = "year", estimator = "ra",
      controls = NULL, ps_controls = NULL, aggregate = "none",
      control_group = "not_yet_treated", graph = FALSE,
      ri = FALSE, rireps = 1000L, ri_method = "bootstrap",
      exclude_pre_periods = 0L, Q = 4L, mode = "common_timing"
    ),
    class = "lwdid_data"
  )
})

test_that("Layer 6: vce=NULL + cluster_var → warning", {
  expect_warning(
    .validate_cross_param_consistency(
      vce = NULL, cluster_var = "cl", rolling = "demean",
      season_var = NULL, tvar = "year", estimator = "ra",
      controls = NULL, ps_controls = NULL, aggregate = "none",
      control_group = "not_yet_treated", graph = FALSE,
      ri = FALSE, rireps = 1000L, ri_method = "bootstrap",
      exclude_pre_periods = 0L, Q = 4L, mode = "common_timing"
    ),
    class = "lwdid_data"
  )
})


# ============================================================================
# Group 6: Layer 7-8 — Data conversion + missing values
# ============================================================================

test_that("Layer 7: original data.frame not modified", {
  df <- data.frame(id = 1:3, year = 2001:2003, y = 1:3)
  df_copy <- data.frame(id = 1:3, year = 2001:2003, y = 1:3)
  dt <- .convert_to_datatable(df)
  data.table::set(dt, j = "y", value = 99)
  expect_equal(df, df_copy)
})

test_that("Layer 7: data.table input also copied", {
  dt_orig <- data.table::data.table(id = 1:3, y = 1:3)
  dt_copy <- .convert_to_datatable(dt_orig)
  data.table::set(dt_copy, j = "y", value = 99)
  expect_equal(dt_orig$y, 1:3)
})

test_that("Layer 7: string ivar converted to numeric", {
  dt <- data.table::data.table(
    id = c("CA", "NY", "TX"),
    y = 1:3
  )
  result <- .convert_string_id(dt, "id")
  expect_true(is.numeric(result$data$id))
  expect_false(is.null(result$id_mapping))
})

test_that("Layer 7: string ID uses alphabetical ordering", {
  dt <- data.table::data.table(
    id = c("TX", "CA", "NY"),
    y = 1:3
  )
  result <- .convert_string_id(dt, "id")
  # CA=1, NY=2, TX=3 (alphabetical)
  expect_equal(result$data$id, c(3, 1, 2))
  expect_equal(result$id_mapping$original_to_numeric$CA, 1)
  expect_equal(result$id_mapping$original_to_numeric$NY, 2)
  expect_equal(result$id_mapping$original_to_numeric$TX, 3)
})

test_that("Layer 7: id_mapping bidirectional correctness", {
  dt <- data.table::data.table(
    id = c("CA", "NY", "TX"),
    y = 1:3
  )
  result <- .convert_string_id(dt, "id")
  o2n <- result$id_mapping$original_to_numeric
  n2o <- result$id_mapping$numeric_to_original
  for (nm in names(o2n)) {
    expect_equal(n2o[[as.character(o2n[[nm]])]], nm)
  }
})

test_that("Layer 7: numeric ivar not converted", {
  dt <- data.table::data.table(id = 1:3, y = 1:3)
  result <- .convert_string_id(dt, "id")
  expect_null(result$id_mapping)
  expect_equal(result$data$id, 1:3)
})

test_that("Layer 7: factor ivar converted", {
  dt <- data.table::data.table(
    id = factor(c("CA", "NY")),
    y = 1:2
  )
  result <- .convert_string_id(dt, "id")
  expect_true(is.numeric(result$data$id))
  expect_false(is.null(result$id_mapping))
})


test_that("Layer 8: NA rows dropped with warning", {
  dt <- data.table::data.table(
    id = c(1, 1, 2, 2),
    year = c(2001, 2002, 2001, 2002),
    y = c(1, NA, 3, 4)
  )
  expect_warning(
    result <- .handle_missing_values(
      dt, "y", "id", "year", "staggered", NULL, NULL
    ),
    class = "lwdid_data"
  )
  expect_equal(nrow(result), 3L)
})

test_that("Layer 8: no NA → no warning", {
  dt <- data.table::data.table(
    id = c(1, 1, 2, 2),
    year = c(2001, 2002, 2001, 2002),
    y = c(1, 2, 3, 4)
  )
  expect_silent(
    result <- .handle_missing_values(
      dt, "y", "id", "year", "staggered", NULL, NULL
    )
  )
  expect_equal(nrow(result), 4L)
})

test_that("Layer 8: CT mode checks d_ and post_ NA", {
  dt <- data.table::data.table(
    id = c(1, 1, 2, 2),
    year = c(2001, 2002, 2001, 2002),
    y = c(1, 2, 3, 4),
    d_ = c(1L, 1L, NA_integer_, 0L),
    post_ = c(0L, 1L, 0L, 1L)
  )
  expect_warning(
    result <- .handle_missing_values(
      dt, "y", "id", "year", "common_timing", "d", "post"
    ),
    class = "lwdid_data"
  )
  expect_equal(nrow(result), 3L)
})


# ============================================================================
# Group 7: Layer 9 — Time index + cohort extraction
# ============================================================================

test_that("Layer 9: annual tindex from 1", {
  dt <- data.table::data.table(
    year = c(2001, 2002, 2003, 2001, 2002, 2003),
    id = c(1, 1, 1, 2, 2, 2)
  )
  result <- .create_time_index(dt, "year")
  expect_equal(sort(unique(result$data$tindex)), 1:3)
  expect_false(result$is_quarterly)
})

test_that("Layer 9: quarterly tq formula", {
  dt <- data.table::data.table(
    year = c(2000, 2000, 2000, 2000),
    quarter = c(1, 2, 3, 4),
    id = rep(1, 4)
  )
  result <- .create_time_index(dt, c("year", "quarter"))
  # tq = (2000-1960)*4 + quarter = 160 + 1..4 = 161..164
  expect_equal(result$data$tq, c(161, 162, 163, 164))
  expect_equal(result$data$tindex, 1:4)
  expect_true(result$is_quarterly)
})

test_that("Layer 9: tindex after dropna uses cleaned min", {
  # If row with min year has NA in y, it gets dropped
  # tindex should still start from 1 based on remaining data
  dt <- data.table::data.table(
    id = c(1, 1, 1, 2, 2, 2),
    year = c(2000, 2001, 2002, 2000, 2001, 2002),
    y = c(NA, 2, 3, 4, 5, 6)
  )
  # Drop NA rows first
  dt_clean <- dt[!is.na(y)]
  result <- .create_time_index(dt_clean, "year")
  # min year in cleaned data is 2000 (unit 2 still has it)
  expect_equal(min(result$data$tindex), 1L)
})


test_that("Layer 9: staggered gvar non-numeric → error", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3),
    year = rep(2001:2003, 2),
    gvar = rep(c("a", "b"), each = 3),
    y = rnorm(6),
    tindex = rep(1:3, 2)
  )
  expect_error(
    .validate_staggered_data(dt, "gvar", "id", "year", "y"),
    class = "lwdid_invalid_staggered_data"
  )
})

test_that("Layer 9: staggered gvar time-varying → error", {
  dt <- data.table::data.table(
    id = c(1, 1, 1),
    year = c(2001, 2002, 2003),
    gvar = c(2002, 2003, 2002),  # varies within unit
    y = rnorm(3),
    tindex = 1:3
  )
  expect_error(
    .validate_staggered_data(dt, "gvar", "id", "year", "y"),
    class = "lwdid_invalid_staggered_data"
  )
})

test_that("Layer 9: staggered gvar negative → error", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3),
    year = rep(2001:2003, 2),
    gvar = rep(c(-1, 2002), each = 3),
    y = rnorm(6),
    tindex = rep(1:3, 2)
  )
  expect_error(
    .validate_staggered_data(dt, "gvar", "id", "year", "y"),
    class = "lwdid_invalid_staggered_data"
  )
})

test_that("Layer 9: staggered gvar -Inf → error", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3),
    year = rep(2001:2003, 2),
    gvar = rep(c(-Inf, 2002), each = 3),
    y = rnorm(6),
    tindex = rep(1:3, 2)
  )
  expect_error(
    .validate_staggered_data(dt, "gvar", "id", "year", "y"),
    class = "lwdid_invalid_staggered_data"
  )
})

test_that("Layer 9: all units NT → error", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    year = rep(2001:2003, 3),
    gvar = rep(c(0, Inf, NA), each = 3),
    y = rnorm(9),
    tindex = rep(1:3, 3)
  )
  expect_error(
    .validate_staggered_data(dt, "gvar", "id", "year", "y"),
    class = "lwdid_invalid_staggered_data"
  )
})

test_that("Layer 9: NT identification", {
  # NA → NT, 0 → NT, +Inf → NT, near-zero → NT
  expect_true(is_never_treated(NA))
  expect_true(is_never_treated(NaN))
  expect_true(is_never_treated(0))
  expect_true(is_never_treated(Inf))
  expect_true(is_never_treated(1e-11))
  expect_false(is_never_treated(2005))
  expect_false(is_never_treated(1))
})

test_that("Layer 9: staggered cohorts sorted correctly", {
  dt <- data.table::data.table(
    id = rep(1:4, each = 5),
    year = rep(2001:2005, 4),
    gvar = rep(c(0, 2004, 2003, 2005), each = 5),
    y = rnorm(20),
    tindex = rep(1:5, 4)
  )
  result <- .validate_staggered_data(dt, "gvar", "id", "year", "y")
  expect_equal(result$cohorts, c(2003L, 2004L, 2005L))
  expect_equal(result$n_cohorts, 3L)
  expect_equal(result$n_never_treated, 1L)
  expect_true(result$has_never_treated)
})


# ============================================================================
# Group 8: Layer 10 — Panel structure
# ============================================================================

test_that("Layer 10: duplicates detected", {
  dt <- data.table::data.table(
    id = c(1, 1, 1, 2, 2),
    year = c(2001, 2001, 2002, 2001, 2002),  # dup for unit 1
    tindex = c(1, 1, 2, 1, 2)
  )
  expect_error(
    .validate_no_duplicates(dt, "id", "year"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 10: no duplicates passes", {
  dt <- data.table::data.table(
    id = c(1, 1, 2, 2),
    tindex = c(1, 2, 1, 2)
  )
  expect_silent(.validate_no_duplicates(dt, "id", "year"))
})

test_that("Layer 10: CT no pre-period → error", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 2),
    tindex = rep(1:2, 3),
    d_ = rep(c(1, 1, 0), each = 2),
    post_ = rep(1L, 6)  # all post
  )
  expect_error(
    .validate_common_timing_structure(dt, "id"),
    class = "lwdid_insufficient_data"
  )
})

test_that("Layer 10: CT no post-period → error", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 2),
    tindex = rep(1:2, 3),
    d_ = rep(c(1, 1, 0), each = 2),
    post_ = rep(0L, 6)  # all pre
  )
  expect_error(
    .validate_common_timing_structure(dt, "id"),
    class = "lwdid_insufficient_data"
  )
})

test_that("Layer 10: CT post varies across units → error", {
  dt <- data.table::data.table(
    id = c(1, 1, 2, 2),
    tindex = c(1, 2, 1, 2),
    d_ = c(1L, 1L, 0L, 0L),
    post_ = c(0L, 1L, 0L, 0L)  # unit 2 has post=0 at t=2
  )
  expect_error(
    .validate_common_timing_structure(dt, "id"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Layer 10: CT post non-monotone → error", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 4),
    tindex = rep(1:4, 2),
    d_ = rep(c(1L, 0L), each = 4),
    post_ = rep(c(0L, 1L, 0L, 1L), 2)  # reversal at t=3
  )
  # post varies at t=3 (0 after 1) — this will be caught by
  # cross-unit consistency first since all units have same pattern
  # Actually post_by_time consistency passes (all units same),
  # but monotonicity fails
  expect_error(
    .validate_common_timing_structure(dt, "id"),
    class = "lwdid_time_discontinuity"
  )
})

test_that("Layer 10: CT no treated units → error", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    tindex = rep(1:4, 3),
    d_ = rep(0L, 12),  # all control
    post_ = rep(c(0L, 0L, 1L, 1L), 3)
  )
  expect_error(
    .validate_common_timing_structure(dt, "id"),
    class = "lwdid_no_treated"
  )
})

test_that("Layer 10: CT no control units → error", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    tindex = rep(1:4, 3),
    d_ = rep(1L, 12),  # all treated
    post_ = rep(c(0L, 0L, 1L, 1L), 3)
  )
  expect_error(
    .validate_common_timing_structure(dt, "id"),
    class = "lwdid_no_control"
  )
})

test_that("Layer 10: CT single treated unit → warning", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    tindex = rep(1:4, 3),
    d_ = rep(c(1L, 0L, 0L), each = 4),
    post_ = rep(c(0L, 0L, 1L, 1L), 3)
  )
  expect_warning(
    .validate_common_timing_structure(dt, "id"),
    class = "lwdid_small_sample"
  )
})

test_that("Layer 10: CT single control unit → warning", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    tindex = rep(1:4, 3),
    d_ = rep(c(1L, 1L, 0L), each = 4),
    post_ = rep(c(0L, 0L, 1L, 1L), 3)
  )
  expect_warning(
    .validate_common_timing_structure(dt, "id"),
    class = "lwdid_small_sample"
  )
})

test_that("Layer 10: CT K and tpost1 correct", {
  dt <- data.table::data.table(
    id = rep(1:4, each = 6),
    tindex = rep(1:6, 4),
    d_ = rep(c(1L, 1L, 0L, 0L), each = 6),
    post_ = rep(c(0L, 0L, 0L, 1L, 1L, 1L), 4)
  )
  result <- .validate_common_timing_structure(dt, "id")
  expect_equal(result$K, 3L)
  expect_equal(result$tpost1, 4L)
  expect_equal(result$n_treated, 2L)
  expect_equal(result$n_control, 2L)
})


test_that("Layer 10: panel balanced → TRUE", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    year = rep(2001:2004, 3)
  )
  expect_true(
    .validate_panel_balance(dt, "id", "warn", "common_timing",
                            NULL, "year")
  )
})

test_that("Layer 10: unbalanced + error → lwdid_unbalanced_panel", {
  dt <- data.table::data.table(
    id = c(rep(1, 4), rep(2, 3), rep(3, 4)),
    year = c(2001:2004, 2001:2003, 2001:2004)
  )
  expect_error(
    .validate_panel_balance(dt, "id", "error", "common_timing",
                            NULL, "year"),
    class = "lwdid_unbalanced_panel"
  )
})

test_that("Layer 10: unbalanced + warn → lwdid_data warning", {
  dt <- data.table::data.table(
    id = c(rep(1, 4), rep(2, 3), rep(3, 4)),
    year = c(2001:2004, 2001:2003, 2001:2004)
  )
  expect_warning(
    .validate_panel_balance(dt, "id", "warn", "common_timing",
                            NULL, "year"),
    class = "lwdid_data"
  )
})

test_that("Layer 10: unbalanced + ignore → silent", {
  dt <- data.table::data.table(
    id = c(rep(1, 4), rep(2, 3), rep(3, 4)),
    year = c(2001:2004, 2001:2003, 2001:2004)
  )
  expect_silent(
    .validate_panel_balance(dt, "id", "ignore", "common_timing",
                            NULL, "year")
  )
})

# ============================================================================
# Group 8a: Deferred cross-parameter checks
# ============================================================================

test_that("Deferred: aggregate=cohort + no NT → error", {
  stag_info <- list(has_never_treated = FALSE, n_never_treated = 0L)
  expect_error(
    .validate_deferred_staggered_checks(
      stag_info, "cohort", "not_yet_treated", NULL,
      data.table::data.table(), "id", "gvar"
    ),
    class = "lwdid_no_never_treated"
  )
})

test_that("Deferred: aggregate=overall + no NT → error", {
  stag_info <- list(has_never_treated = FALSE, n_never_treated = 0L)
  expect_error(
    .validate_deferred_staggered_checks(
      stag_info, "overall", "not_yet_treated", NULL,
      data.table::data.table(), "id", "gvar"
    ),
    class = "lwdid_no_never_treated"
  )
})

test_that("Deferred: control_group=never_treated + no NT → error", {
  stag_info <- list(has_never_treated = FALSE, n_never_treated = 0L)
  expect_error(
    .validate_deferred_staggered_checks(
      stag_info, "none", "never_treated", NULL,
      data.table::data.table(), "id", "gvar"
    ),
    class = "lwdid_no_never_treated"
  )
})

test_that("Deferred: aggregate + not_yet_treated → auto-switch + warning", {
  stag_info <- list(has_never_treated = TRUE, n_never_treated = 3L)
  dt <- data.table::data.table(
    id = rep(1:4, each = 3),
    gvar = rep(c(0, 2003, 2004, 2005), each = 3)
  )
  expect_warning(
    result <- .validate_deferred_staggered_checks(
      stag_info, "cohort", "not_yet_treated", NULL,
      dt, "id", "gvar"
    ),
    class = "lwdid_control_group_switch"
  )
  expect_equal(result, "never_treated")
})

test_that("Deferred: control_group=auto + has NT → never_treated", {
  stag_info <- list(has_never_treated = TRUE, n_never_treated = 2L)
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    gvar = rep(c(0, 2003, 2004), each = 3)
  )
  result <- .validate_deferred_staggered_checks(
    stag_info, "none", "auto", NULL, dt, "id", "gvar"
  )
  expect_equal(result, "never_treated")
})

test_that("Deferred: control_group=auto + no NT → not_yet_treated", {
  stag_info <- list(has_never_treated = FALSE, n_never_treated = 0L)
  dt <- data.table::data.table(
    id = rep(1:2, each = 3),
    gvar = rep(c(2003, 2004), each = 3)
  )
  result <- .validate_deferred_staggered_checks(
    stag_info, "none", "auto", NULL, dt, "id", "gvar"
  )
  expect_equal(result, "not_yet_treated")
})


# ============================================================================
# Group 5a: gid validation
# ============================================================================

test_that("gid=NULL passes (skips validation)", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    d_ = rep(c(1L, 1L, 0L), each = 4)
  )
  expect_silent(.validate_gid_treated(NULL, dt, "id", "d_"))
})

test_that("gid as treated unit passes", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    d_ = rep(c(1L, 1L, 0L), each = 4)
  )
  expect_silent(.validate_gid_treated(1, dt, "id", "d_"))
})

test_that("gid as control unit → error", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    d_ = rep(c(1L, 1L, 0L), each = 4)
  )
  expect_error(
    .validate_gid_treated(3, dt, "id", "d_"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("gid not in data → error", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    d_ = rep(c(1L, 1L, 0L), each = 4)
  )
  expect_error(
    .validate_gid_treated(99, dt, "id", "d_"),
    class = "lwdid_invalid_parameter"
  )
})

# ============================================================================
# Group 9: Pre-period sufficiency
# ============================================================================

test_that("Pre-period: CT demean + 1 pre → passes", {
  expect_silent(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "demean",
      n_pre_periods = 1L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    )
  )
})

test_that("Pre-period: CT demean + 0 pre → error", {
  expect_error(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "demean",
      n_pre_periods = 0L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    ),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("Pre-period: CT detrend + 1 pre → error", {
  expect_error(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "detrend",
      n_pre_periods = 1L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    ),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("Pre-period: CT detrend + 2 pre → passes", {
  expect_silent(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "detrend",
      n_pre_periods = 2L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    )
  )
})

test_that("Pre-period: CT demeanq needs Q+1", {
  # Q=4, need 5 pre-periods
  expect_error(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "demeanq",
      n_pre_periods = 4L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    ),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_silent(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "demeanq",
      n_pre_periods = 5L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    )
  )
})

test_that("Pre-period: CT detrendq needs Q+2", {
  # Q=4, need 6 pre-periods
  expect_error(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "detrendq",
      n_pre_periods = 5L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    ),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_silent(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "detrendq",
      n_pre_periods = 6L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    )
  )
})

test_that("Pre-period: exclude_pre_periods reduces available", {
  # detrend needs 2, have 3, exclude 2 → 1 left → error
  expect_error(
    .validate_pre_period_sufficiency(
      mode = "common_timing", rolling = "detrend",
      n_pre_periods = 3L, exclude_pre_periods = 2L,
      Q = 4L, cohorts = NULL, T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    ),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("Pre-period: staggered all cohorts insufficient → error", {
  expect_error(
    .validate_pre_period_sufficiency(
      mode = "staggered", rolling = "detrend",
      n_pre_periods = 0L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = c(2002L, 2003L), T_min = 2002L,
      dt = NULL, ivar = NULL, gvar = NULL
    ),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("Pre-period: staggered some cohorts insufficient → warning", {
  # cohort 2002: 2002-2001=1 pre (ok for demean)
  # cohort 2001: 2001-2001=0 pre (insufficient)
  expect_warning(
    .validate_pre_period_sufficiency(
      mode = "staggered", rolling = "demean",
      n_pre_periods = 0L, exclude_pre_periods = 0L,
      Q = 4L, cohorts = c(2001L, 2002L), T_min = 2001L,
      dt = NULL, ivar = NULL, gvar = NULL
    ),
    class = "lwdid_data"
  )
})

# ============================================================================
# Group 5a (cont.): Staggered gid validation via deferred checks
# ============================================================================
# NOTE: Column name "first_treat" used instead of "gvar" to avoid
# data.table scoping conflict with the function parameter `gvar`.

test_that("Staggered gid: gid not found in data → error", {
  stag_info <- list(has_never_treated = TRUE, n_never_treated = 1L)
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    year = rep(2001:2003, 3),
    first_treat = rep(c(0, 2002, 2003), each = 3)
  )
  expect_error(
    .validate_deferred_staggered_checks(
      stag_info, "none", "auto", gid = 99,
      dt, "id", "first_treat"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Staggered gid: gid is never-treated unit (gvar=0) → error", {
  # Unit 1 has first_treat=0 → never-treated
  stag_info <- list(has_never_treated = TRUE, n_never_treated = 1L)
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    year = rep(2001:2003, 3),
    first_treat = rep(c(0, 2002, 2003), each = 3)
  )
  expect_error(
    .validate_deferred_staggered_checks(
      stag_info, "none", "auto", gid = 1,
      dt, "id", "first_treat"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Staggered gid: gid is valid treated unit → passes", {
  stag_info <- list(has_never_treated = TRUE, n_never_treated = 1L)
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    year = rep(2001:2003, 3),
    first_treat = rep(c(0, 2002, 2003), each = 3)
  )
  # gid=2 is treated (first_treat=2002), should pass without error
  result <- .validate_deferred_staggered_checks(
    stag_info, "none", "auto", gid = 2,
    dt, "id", "first_treat"
  )
  expect_true(result %in% c("never_treated", "not_yet_treated"))
})

test_that("Staggered gid: gid with NA gvar (never-treated) → error", {
  stag_info <- list(has_never_treated = TRUE, n_never_treated = 1L)
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    year = rep(2001:2003, 3),
    first_treat = rep(c(NA, 2002, 2003), each = 3)
  )
  expect_error(
    .validate_deferred_staggered_checks(
      stag_info, "none", "auto", gid = 1,
      dt, "id", "first_treat"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("Staggered gid: gid with Inf gvar (never-treated) → error", {
  stag_info <- list(has_never_treated = TRUE, n_never_treated = 1L)
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    year = rep(2001:2003, 3),
    first_treat = rep(c(Inf, 2002, 2003), each = 3)
  )
  expect_error(
    .validate_deferred_staggered_checks(
      stag_info, "none", "auto", gid = 1,
      dt, "id", "first_treat"
    ),
    class = "lwdid_invalid_parameter"
  )
})

# ============================================================================
# Group 10: Season diversity tests (.validate_season_diversity / .check_unit_season_diversity)
# ============================================================================

test_that("Season diversity CT: pre-period < 2 seasons → error", {
  # Unit 1: pre has only quarter=1, post has quarter=2
  dt <- data.table::data.table(
    id = rep(1L, 4),
    quarter = c(1L, 1L, 2L, 2L),
    post_ = c(0L, 0L, 1L, 1L)
  )
  expect_error(
    .validate_season_diversity(dt, "id", "quarter", "post_", Q = 4L),
    class = "lwdid_insufficient_quarter_diversity"
  )
})

test_that("Season diversity CT: post season not covered by pre → error", {

  # Unit 1: pre has Q1,Q2; post has Q1,Q3 → Q3 uncovered
  dt <- data.table::data.table(
    id = rep(1L, 6),
    quarter = c(1L, 2L, 1L, 2L, 1L, 3L),
    post_ = c(0L, 0L, 0L, 0L, 1L, 1L)
  )
  expect_error(
    .validate_season_diversity(dt, "id", "quarter", "post_", Q = 4L),
    class = "lwdid_insufficient_quarter_diversity"
  )
})

test_that("Season diversity CT: all post seasons covered → pass", {
  # Unit 1: pre has Q1,Q2,Q3; post has Q1,Q2
  dt <- data.table::data.table(
    id = rep(1L, 5),
    quarter = c(1L, 2L, 3L, 1L, 2L),
    post_ = c(0L, 0L, 0L, 1L, 1L)
  )
  expect_silent(
    .validate_season_diversity(dt, "id", "quarter", "post_", Q = 4L)
  )
})

test_that("Season diversity Staggered: treated unit pre < 2 seasons → error", {
  # Unit 1 treated (g=2003): times 2001,2002 → pre = 2001,2002 (both same quarter=1)
  dt <- data.table::data.table(
    id = c(1L, 1L, 1L, 1L),
    year = c(2001L, 2002L, 2003L, 2004L),
    quarter = c(1L, 1L, 2L, 2L),
    gvar = c(2003L, 2003L, 2003L, 2003L)
  )
  expect_error(
    .validate_season_diversity(dt, "id", "quarter", NULL, Q = 4L,
                               gvar = "gvar", tvar = "year"),
    class = "lwdid_insufficient_quarter_diversity"
  )
})

test_that("Season diversity Staggered: NT units are skipped", {
  # Unit 1 is NT (g=Inf): only 1 season in pre, but should not error
  # Unit 2 is treated (g=2003): has 2 seasons in pre, post covered
  dt <- data.table::data.table(
    id = c(rep(1L, 4), rep(2L, 4)),
    year = rep(c(2001L, 2002L, 2003L, 2004L), 2),
    quarter = c(1L, 1L, 1L, 1L,   # unit 1: all Q1 (NT, skipped)
                1L, 2L, 1L, 2L),   # unit 2: Q1,Q2 pre; Q1,Q2 post
    gvar = c(rep(Inf, 4), rep(2003L, 4))
  )
  expect_silent(
    .validate_season_diversity(dt, "id", "quarter", NULL, Q = 4L,
                               gvar = "gvar", tvar = "year")
  )
})

test_that("Season diversity: freq_label uses FREQ_LABELS (Q=4 → 'quarter')", {
  dt <- data.table::data.table(
    id = rep(1L, 2),
    quarter = c(1L, 1L),
    post_ = c(0L, 1L)
  )
  err <- tryCatch(
    .validate_season_diversity(dt, "id", "quarter", "post_", Q = 4L),
    lwdid_insufficient_quarter_diversity = function(e) e
  )
  expect_true(grepl("quarter", err$message))
})

test_that("Season diversity: Q=12 uses 'month' label", {
  dt <- data.table::data.table(
    id = rep(1L, 2),
    month = c(1L, 1L),
    post_ = c(0L, 1L)
  )
  err <- tryCatch(
    .validate_season_diversity(dt, "id", "month", "post_", Q = 12L),
    lwdid_insufficient_quarter_diversity = function(e) e
  )
  expect_true(grepl("month", err$message))
})

test_that("Season diversity: Q=7 (non-standard) uses 'season' label", {
  dt <- data.table::data.table(
    id = rep(1L, 2),
    s = c(1L, 1L),
    post_ = c(0L, 1L)
  )
  err <- tryCatch(
    .validate_season_diversity(dt, "id", "s", "post_", Q = 7L),
    lwdid_insufficient_quarter_diversity = function(e) e
  )
  expect_true(grepl("season", err$message))
})

test_that("check_unit_season_diversity: error metadata correct", {
  # Pre has Q1 only → error with unit, missing_quarters, pre_quarters
  dt_pre <- data.table::data.table(quarter = 1L)
  dt_post <- data.table::data.table(quarter = 2L)
  err <- tryCatch(
    .check_unit_season_diversity(42L, dt_pre, dt_post, "quarter", "quarter"),
    lwdid_insufficient_quarter_diversity = function(e) e
  )
  expect_equal(err$unit, 42L)
  expect_equal(err$pre_quarters, 1L)
})

test_that("check_unit_season_diversity: uncovered post seasons in metadata", {
  dt_pre <- data.table::data.table(quarter = c(1L, 2L))
  dt_post <- data.table::data.table(quarter = c(1L, 3L))
  err <- tryCatch(
    .check_unit_season_diversity(1L, dt_pre, dt_post, "quarter", "quarter"),
    lwdid_insufficient_quarter_diversity = function(e) e
  )
  expect_equal(err$missing_quarters, 3L)
  expect_equal(err$pre_quarters, c(1L, 2L))
  expect_equal(err$post_quarters, c(1L, 3L))
})

# ============================================================================
# Group 11: Return value structure tests
# ============================================================================

test_that("Return value: CT mode has all required fields", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  expected_fields <- c(
    "mode", "data", "n_obs", "n_units", "n_periods",
    "K", "tpost1", "T_min", "T_max", "is_quarterly",
    "n_pre_periods", "cohort_pre_periods",
    "cohorts", "n_cohorts", "cohort_sizes",
    "has_nt", "n_treated", "n_control", "n_nt",
    "periods", "balanced", "id_mapping",
    "control_group_resolved", "freq_info", "validated_params"
  )
  for (f in expected_fields) {
    expect_true(f %in% names(result),
                info = paste("Missing field:", f))
  }
})

test_that("Return value: CT mode field types correct", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  expect_equal(result$mode, "common_timing")
  expect_true(data.table::is.data.table(result$data))
  expect_true(is.integer(result$n_obs) || is.numeric(result$n_obs))
  expect_true(is.integer(result$n_units) || is.numeric(result$n_units))
  expect_true(is.integer(result$K) || is.numeric(result$K))
  expect_true(is.integer(result$tpost1) || is.numeric(result$tpost1))
  expect_true(is.integer(result$T_min))
  expect_true(is.integer(result$T_max))
  expect_true(is.logical(result$is_quarterly))
  expect_null(result$cohorts)
  expect_null(result$n_cohorts)
  expect_null(result$cohort_sizes)
  expect_null(result$freq_info)  # auto_detect_frequency=FALSE by default
})

test_that("Return value: Staggered mode has all required fields", {
  stag <- make_stag_data()
  result <- suppressWarnings(validate_inputs(
    data = stag, y = "y", gvar = "gvar", ivar = "id", tvar = "year",
    rolling = "demean", aggregate = "none",
    control_group = "not_yet_treated"
  ))
  expect_equal(result$mode, "staggered")
  expect_true(!is.null(result$cohorts))
  expect_true(!is.null(result$n_cohorts))
  expect_true(!is.null(result$cohort_sizes))
  expect_true(is.logical(result$has_nt))
})

test_that("Return value: validated_params contains all input parameters", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "Demean", estimator = "ra",
    alpha = 0.1
  ))
  vp <- result$validated_params
  expect_equal(vp$depvar, "y")
  expect_equal(vp$ivar, "id")
  expect_equal(vp$d, "d")
  expect_equal(vp$post, "post")
  expect_equal(vp$rolling, "demean")  # lowercased
  expect_equal(vp$estimator, "ra")
  expect_equal(vp$alpha, 0.1)
})

test_that("Return value: rolling is lowercased in validated_params", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "DETREND"
  ))
  expect_equal(result$validated_params$rolling, "detrend")
})

test_that("Return value: control_group_resolved for CT", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean",
    control_group = "not_yet_treated"
  ))
  # CT mode: auto resolves to "not_yet_treated"
  expect_equal(result$control_group_resolved, "not_yet_treated")
})

test_that("Return value: Staggered auto control_group with NT → never_treated", {
  stag <- make_stag_data()  # has NT units (gvar=Inf)
  result <- suppressWarnings(validate_inputs(
    data = stag, y = "y", gvar = "gvar", ivar = "id", tvar = "year",
    rolling = "demean", aggregate = "none",
    control_group = "auto"
  ))
  expect_equal(result$control_group_resolved, "never_treated")
})

test_that("Return value: Staggered auto control_group without NT → not_yet_treated", {
  # Create staggered data with no NT units
  stag_no_nt <- data.frame(
    id = rep(1:6, each = 5),
    year = rep(2001:2005, 6),
    y = rnorm(30),
    gvar = rep(c(2003L, 2003L, 2004L, 2004L, 2005L, 2005L), each = 5)
  )
  result <- suppressWarnings(validate_inputs(
    data = stag_no_nt, y = "y", gvar = "gvar", ivar = "id", tvar = "year",
    rolling = "demean", aggregate = "none",
    control_group = "auto"
  ))
  expect_equal(result$control_group_resolved, "not_yet_treated")
})

# ============================================================================
# Group 12: End-to-end integration tests
# ============================================================================

test_that("E2E: valid CT data passes validate_inputs() completely", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  expect_equal(result$mode, "common_timing")
  expect_true(result$n_obs > 0)
  expect_true(result$n_units >= 3)
  expect_true(result$K > 0)
  expect_true(result$tpost1 > 0)
  expect_true(result$n_treated > 0)
  expect_true(result$n_control > 0)
  # data should be a data.table with tindex column
  expect_true("tindex" %in% names(result$data))
  # d_ and post_ columns should exist
  expect_true("d_" %in% names(result$data))
  expect_true("post_" %in% names(result$data))
})

test_that("E2E: valid Staggered data passes validate_inputs() completely", {
  stag <- make_stag_data()
  result <- suppressWarnings(validate_inputs(
    data = stag, y = "y", gvar = "gvar", ivar = "id", tvar = "year",
    rolling = "demean", aggregate = "none",
    control_group = "not_yet_treated"
  ))
  expect_equal(result$mode, "staggered")
  expect_true(result$n_obs > 0)
  expect_true(result$n_units >= 3)
  expect_true(!is.null(result$cohorts))
  expect_true(length(result$cohorts) > 0)
  expect_true("tindex" %in% names(result$data))
  # No d_/post_ columns in staggered mode
  expect_false("d_" %in% names(result$data))
  expect_false("post_" %in% names(result$data))
})

test_that("E2E: quarterly CT data passes validate_inputs()", {
  # Create quarterly panel: tvar = c("year", "quarter")
  n_units <- 6
  years <- 2000:2003
  quarters <- 1:4
  grid <- expand.grid(quarter = quarters, year = years, id = 1:n_units)
  grid$y <- rnorm(nrow(grid))
  # Treatment: units 1-3 treated, 4-6 control
  grid$d <- ifelse(grid$id <= 3, 1L, 0L)
  # Post: year >= 2002
  grid$post <- ifelse(grid$year >= 2002, 1L, 0L)
  result <- suppressWarnings(validate_inputs(
    data = grid, y = "y", d = "d", ivar = "id",
    tvar = c("year", "quarter"), post = "post", rolling = "demean"
  ))
  expect_equal(result$mode, "common_timing")
  expect_true(result$is_quarterly)
  expect_true("tq" %in% names(result$data))
  expect_true("tindex" %in% names(result$data))
  # Verify tq calculation: (year - 1960) * 4 + quarter
  first_row <- result$data[1]
  expected_tq <- (first_row$year - 1960L) * 4L + first_row$quarter
  expect_equal(first_row$tq, expected_tq)
})

test_that("E2E: CT numerical values are reasonable", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  # make_ct_data defaults: n_units=10, n_periods=8, K=4
  # years 2001-2008, post after year 2004 (i.e., 2005-2008)
  # K = 4 (pre-periods: 2001-2004), tpost1 = 5
  expect_equal(result$n_units, 10L)
  expect_equal(result$n_periods, 8L)
  expect_equal(result$n_obs, 80L)  # 10 units * 8 years
  expect_equal(result$K, 4L)
  expect_equal(result$tpost1, 5L)
  expect_equal(result$T_min, 2001L)
  expect_equal(result$T_max, 2008L)
  expect_equal(result$n_treated, 5L)  # units 1-5
  expect_equal(result$n_control, 5L)  # units 6-10
})

test_that("E2E: Staggered numerical values are reasonable", {
  stag <- make_stag_data()
  result <- suppressWarnings(validate_inputs(
    data = stag, y = "y", gvar = "gvar", ivar = "id", tvar = "year",
    rolling = "demean", aggregate = "none",
    control_group = "not_yet_treated"
  ))
  # make_stag_data defaults: n_units=12, n_periods=10, cohorts=c(2005,2007), n_nt=4
  # years 2001-2010, units 1-4 NT (gvar=0), 5-8 cohort 2005, 9-12 cohort 2007
  expect_equal(result$n_units, 12L)
  expect_equal(result$T_min, 2001L)
  expect_equal(result$T_max, 2010L)
  expect_equal(sort(result$cohorts), c(2005L, 2007L))
  expect_equal(result$n_cohorts, 2L)
  expect_true(result$has_nt)
  expect_equal(result$n_nt, 4L)  # units 1-4
  expect_equal(result$n_treated, 8L)  # units 5-12
})

test_that("E2E: error condition class names are consistent with conditions.R", {
  # Verify that errors thrown use classes defined in conditions.R
  # Test a few representative errors and check their class hierarchy

  # lwdid_invalid_parameter
  err1 <- tryCatch(
    validate_inputs(
      data = data.frame(), y = "y", d = "d", ivar = "id",
      tvar = "year", post = "post"
    ),
    lwdid_invalid_parameter = function(e) e
  )
  expect_true(inherits(err1, "lwdid_invalid_parameter"))
  expect_true(inherits(err1, "lwdid_error"))
  expect_true(inherits(err1, "error"))

  # lwdid_missing_column
  ct <- make_ct_data()
  err2 <- tryCatch(
    validate_inputs(
      data = ct, y = "nonexistent", d = "d", ivar = "id",
      tvar = "year", post = "post"
    ),
    lwdid_missing_column = function(e) e
  )
  expect_true(inherits(err2, "lwdid_missing_column"))
  expect_true(inherits(err2, "lwdid_error"))
})

# ============================================================================
# Group 13: Frequency detection tests (.detect_frequency / .detect_frequency_numeric)
# ============================================================================

test_that("Frequency detection: annual data (1 obs/year) → annual, Q=1", {
  # 5 units, 1 obs per year, years 2001-2010
  dt <- data.table::data.table(
    id = rep(1:5, each = 10),
    year = rep(2001:2010, 5),
    y = rnorm(50)
  )
  result <- .detect_frequency(dt, "year", "id")
  expect_equal(result$frequency, "annual")
  expect_equal(result$Q, 1L)
  expect_true(result$confidence > 0.5)
  expect_equal(result$method, "obs_per_year")
})

test_that("Frequency detection: quarterly data (4 obs/year) → quarterly, Q=4", {
  # 3 units, 4 obs per year, years 2001-2005
  dt <- data.table::data.table(
    id = rep(1:3, each = 20),
    year = rep(rep(2001:2005, each = 4), 3),
    quarter = rep(rep(1:4, 5), 3),
    y = rnorm(60)
  )
  # Use year as tvar (4 obs per year value)
  result <- .detect_frequency(dt, "year", "id")
  expect_equal(result$frequency, "quarterly")
  expect_equal(result$Q, 4L)
  expect_true(result$confidence > 0.5)
})

test_that("Frequency detection: monthly data (12 obs/year) → monthly, Q=12", {
  # 2 units, 12 obs per year, years 2001-2003
  dt <- data.table::data.table(
    id = rep(1:2, each = 36),
    year = rep(rep(2001:2003, each = 12), 2),
    month = rep(rep(1:12, 3), 2),
    y = rnorm(72)
  )
  result <- .detect_frequency(dt, "year", "id")
  expect_equal(result$frequency, "monthly")
  expect_equal(result$Q, 12L)
  expect_true(result$confidence > 0.5)
})

test_that("Frequency detection: insufficient data (< 2 time values) → warning", {
  dt <- data.table::data.table(
    year = 2001L,
    y = 1.0
  )
  expect_warning(
    result <- .detect_frequency(dt, "year"),
    class = "lwdid_data"
  )
  expect_null(result$frequency)
  expect_null(result$Q)
})

test_that("Frequency detection: missing tvar column → warning", {
  dt <- data.table::data.table(
    id = 1:5,
    y = rnorm(5)
  )
  expect_warning(
    result <- .detect_frequency(dt, "nonexistent"),
    class = "lwdid_data"
  )
  expect_null(result$frequency)
})

test_that("Frequency detection: consecutive integer index → low confidence + warning", {
  # Non-year-like values (1-20), interval=1 → ambiguous
  dt <- data.table::data.table(
    id = rep(1:2, each = 10),
    t = rep(1:10, 2),
    y = rnorm(20)
  )
  expect_warning(
    result <- .detect_frequency(dt, "t", "id"),
    class = "lwdid_data"
  )
  expect_equal(result$confidence, 0.3)
})

test_that("Frequency detection: auto_detect_frequency=FALSE skips detection in validate_inputs", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean",
    auto_detect_frequency = FALSE
  ))
  expect_null(result$freq_info)
})

test_that("Frequency detection: auto_detect_frequency=TRUE returns freq_info", {
  ct <- make_ct_data()
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean",
    auto_detect_frequency = TRUE
  ))
  expect_true(!is.null(result$freq_info))
  expect_equal(result$freq_info$frequency, "annual")
  expect_equal(result$freq_info$Q, 1L)
})

test_that("Frequency detection: without ivar uses global analysis", {
  dt <- data.table::data.table(
    year = rep(2001:2010, 1),
    y = rnorm(10)
  )
  result <- .detect_frequency(dt, "year")
  expect_equal(result$frequency, "annual")
  expect_equal(result$Q, 1L)
})

test_that("Frequency detection: details contain t_min, t_max, t_range", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 5),
    year = rep(2001:2005, 2),
    y = rnorm(10)
  )
  result <- .detect_frequency(dt, "year", "id")
  expect_equal(result$details$t_min, 2001)
  expect_equal(result$details$t_max, 2005)
  expect_equal(result$details$t_range, 4)
})

# ============================================================================
# Group 14: Cohort float comparison tests (.get_cohort_mask)
# ============================================================================

test_that("get_cohort_mask: exact integer match → TRUE", {
  mask <- .get_cohort_mask(c(2003L, 2005L, 2003L, 2007L), 2003L)
  expect_equal(mask, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("get_cohort_mask: floating-point error within tolerance → TRUE", {
  # Simulate floating-point arithmetic: 2003 + tiny error
  gvar_vals <- c(2003 + 1e-10, 2005, 2003 - 1e-10, 2007)
  mask <- .get_cohort_mask(gvar_vals, 2003)
  expect_equal(mask, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("get_cohort_mask: no match → all FALSE", {
  mask <- .get_cohort_mask(c(2003, 2005, 2007), 2009)
  expect_equal(mask, c(FALSE, FALSE, FALSE))
})

test_that("get_cohort_mask: tolerance boundary — diff >= 1e-9 → FALSE", {
  # COHORT_FLOAT_TOLERANCE = 1e-9
  # At scale 2003, 2003+1e-9-2003 ≈ 9.9999e-10 due to floating-point,
  # so we use 1.1e-9 offset which produces a true diff > 1e-9
  mask <- .get_cohort_mask(2003 + 1.1e-9, 2003)
  expect_false(mask)
  # Verify the actual difference is indeed >= 1e-9
  actual_diff <- abs(2003 + 1.1e-9 - 2003)
  expect_true(actual_diff >= 1e-9)
})

test_that("get_cohort_mask: diff just under tolerance → TRUE", {
  # abs(2003 + 9.99e-10 - 2003) = 9.99e-10 < 1e-9
  mask <- .get_cohort_mask(2003 + 9.99e-10, 2003)
  expect_true(mask)
})

test_that("get_cohort_mask: diff clearly above tolerance → FALSE", {
  mask <- .get_cohort_mask(2003 + 1e-8, 2003)
  expect_false(mask)
})

test_that("get_cohort_mask: empty vector → empty logical", {
  mask <- .get_cohort_mask(numeric(0), 2003)
  expect_equal(mask, logical(0))
})

test_that("get_cohort_mask: single element match", {
  mask <- .get_cohort_mask(2005, 2005)
  expect_true(mask)
})

test_that("get_cohort_mask: negative tolerance direction", {
  # 2003 - 5e-10 should match 2003

  mask <- .get_cohort_mask(2003 - 5e-10, 2003)
  expect_true(mask)
})

# ============================================================================
# Group 15: Execution order verification tests
# ============================================================================

test_that("Execution order: binarization on data.table copy — original unchanged", {
  ct <- make_ct_data()
  original_names <- names(ct)
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  # Original data.frame should NOT have d_ or post_ columns
  expect_false("d_" %in% names(ct))
  expect_false("post_" %in% names(ct))
  expect_false("tindex" %in% names(ct))
  # Column names should be unchanged

  expect_equal(names(ct), original_names)
})

test_that("Execution order: binarization before dropna — NA in d handled correctly", {
  # Create CT data where one row has d=NA
  # The binarization should produce d_=NA, then dropna removes that row
  ct <- make_ct_data()
  ct$d[1] <- NA  # unit 1, first year
  n_original <- nrow(ct)

  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  # The row with d=NA should have been removed
  expect_true(result$n_obs < n_original)
  # No NA in d_ column of result data
  expect_false(any(is.na(result$data$d_)))
})

test_that("Execution order: tindex created after dropna — dropped rows don't affect tindex", {
  # Create CT data where the earliest year row for one unit has NA in y
  # After dropna removes it, tindex should still start from 1 based on
  # the remaining data's min time
  n_units <- 5
  years <- 2001:2010
  ct <- data.frame(
    id = rep(1:n_units, each = length(years)),
    year = rep(years, n_units),
    y = rnorm(n_units * length(years)),
    d = rep(c(1L, 1L, 0L, 0L, 0L), each = length(years)),
    post = rep(ifelse(years >= 2006, 1L, 0L), n_units)
  )
  # Set y=NA for ALL units at year 2001 → after dropna, min year = 2002
  ct$y[ct$year == 2001] <- NA

  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  # tindex should start from 1 (based on min year in cleaned data = 2002)
  expect_equal(min(result$data$tindex), 1L)
  # T_min should be 2002 (not 2001)
  expect_equal(result$T_min, 2002L)
  # Number of periods should be 9 (2002-2010)
  expect_equal(result$n_periods, 9L)
})

test_that("Execution order: binarization values are correct in returned data", {
  # Verify that d_ and post_ in returned data are properly binarized
  ct <- make_ct_data()
  # Set some d values to non-standard (but valid) values
  ct$d <- ifelse(ct$d == 1, 5, 0)  # 5 should become 1 after binarization
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  # d_ should only contain 0 and 1
  expect_true(all(result$data$d_ %in% c(0L, 1L)))
  # The treated units should have d_=1
  expect_true(any(result$data$d_ == 1L))
  expect_true(any(result$data$d_ == 0L))
})

test_that("Execution order: string ID conversion happens before panel checks", {
  # Use string IDs — they should be converted to numeric before panel validation
  ct <- data.frame(
    id = rep(c("CA", "NY", "TX", "FL", "OH"), each = 10),
    year = rep(2001:2010, 5),
    y = rnorm(50),
    d = rep(c(1L, 1L, 0L, 0L, 0L), each = 10),
    post = rep(ifelse(2001:2010 >= 2006, 1L, 0L), 5)
  )
  result <- suppressWarnings(validate_inputs(
    data = ct, y = "y", d = "d", ivar = "id", tvar = "year",
    post = "post", rolling = "demean"
  ))
  # id_mapping should exist (string → numeric conversion happened)
  expect_true(!is.null(result$id_mapping))
  # Numeric IDs should be used in returned data
  expect_true(is.numeric(result$data$id))
  # Alphabetical order: CA=1, FL=2, NY=3, OH=4, TX=5
  expect_equal(result$id_mapping$numeric_to_original[["1"]], "CA")
  expect_equal(result$id_mapping$numeric_to_original[["2"]], "FL")
  expect_equal(result$id_mapping$numeric_to_original[["3"]], "NY")
})
