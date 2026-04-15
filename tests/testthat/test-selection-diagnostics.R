library(testthat)

make_selection_observed_y_panel <- function() {
  data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    gvar = c(2, 2, 0, 0),
    y = c(NA_real_, 5, 3, 4)
  )
}

make_selection_attrition_balance_panel <- function() {
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

make_selection_no_gvar_panel <- function() {
  data.frame(
    id = rep(1:2, each = 3),
    time = rep(1:3, times = 2),
    y = c(1, 2, 3, 4, 5, 6)
  )
}

make_selection_all_missing_panel <- function() {
  data.frame(
    id = rep(1:2, each = 3),
    time = rep(1:3, times = 2),
    gvar = rep(c(2, 0), each = 3),
    y = rep(NA_real_, 6)
  )
}

make_selection_common_timing_panel <- function() {
  data.frame(
    id = c(1, 1, 1, 2, 2, 2),
    time = c(1, 2, 3, 1, 2, 3),
    d = c(0, 0, 1, 0, 0, 0),
    y = c(1, NA_real_, 5, 2, 3, 4)
  )
}

make_selection_single_unit_panel <- function() {
  data.frame(
    id = 1L,
    time = 1:3,
    gvar = 2L,
    y = c(1, 2, 3)
  )
}

make_selection_covariate_mar_panel <- function() {
  data.frame(
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
  )
}

make_selection_layer5_attrition_panel <- function() {
  unit_specs <- data.frame(
    id = 1:12,
    gvar = c(rep(6L, 4), rep(8L, 4), rep(0L, 4)),
    last_obs = c(10L, 8L, 7L, 10L, 10L, 9L, 10L, 10L, 10L, 10L, 10L, 10L),
    stringsAsFactors = FALSE
  )

  panel <- do.call(
    rbind,
    lapply(seq_len(nrow(unit_specs)), function(i) {
      spec <- unit_specs[i, ]
      data.frame(
        id = spec$id,
        time = seq_len(spec$last_obs),
        gvar = spec$gvar,
        y = spec$id * 100 + seq_len(spec$last_obs)
      )
    })
  )

  panel$y[panel$id == 12L & panel$time == 10L] <- NA_real_
  panel
}

test_that("E8-05.1: .is_never_treated honors custom markers and both infinities", {
  expect_equal(
    .is_never_treated(
      c(NA_real_, 0, Inf, -Inf, 99, 5),
      never_treated_values = c(0, 99, Inf)
    ),
    c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  )
})

test_that("E8-05.1: missing-rate helpers use observed Y availability on full grid", {
  panel <- make_selection_observed_y_panel()

  rates <- .compute_missing_rates(panel, ivar = "id", tvar = "time", y = "y")
  by_period <- .compute_missing_by_period(panel, ivar = "id", tvar = "time", y = "y")
  by_cohort <- .compute_missing_by_cohort(
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
  expect_equal(names(by_cohort), "2")
  expect_equal(unname(by_cohort[["2"]]), 0.5, tolerance = 1e-12)
})

test_that("E8-05 observed-Y unit stats exclude y-NA pre-periods from usability", {
  panel <- make_selection_observed_y_panel()

  stats <- .compute_unit_stats(
    panel,
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    y = "y",
    never_treated_values = c(0, Inf)
  )

  treated <- stats[stats$unit_id == 1, ]
  never_treated <- stats[stats$unit_id == 2, ]

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
  expect_equal(treated$reason_if_excluded, "No pre-treatment observations")

  expect_false(never_treated$is_treated)
  expect_true(is.na(never_treated$cohort))
})

test_that("E8-05.3: attrition analysis uses observed Y availability for dropout timing", {
  panel <- make_selection_attrition_balance_panel()

  attrition <- .compute_attrition_analysis(
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

test_that("E8-05.3: balance statistics use observed Y pre-period usability", {
  panel <- make_selection_attrition_balance_panel()

  balance <- .compute_balance_statistics(
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

test_that("E8-05.6: diagnose_selection_mechanism works without gvar in common-timing mode", {
  panel <- make_selection_no_gvar_panel()

  diagnosis <- diagnose_selection_mechanism(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = NULL,
    verbose = FALSE
  )
  unit_stats <- get_unit_missing_stats(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = NULL
  )

  expect_s3_class(diagnosis, "lwdid_selection_diagnosis")
  expect_true(diagnosis$balance_statistics$is_balanced)
  expect_identical(diagnosis$missing_rate_by_cohort, list())
  expect_identical(toupper(diagnosis$selection_risk), "LOW")
  expect_equal(diagnosis$selection_risk_score, 0)
  expect_true(all(is.na(unit_stats$cohort)))
  expect_false(any(unit_stats$is_treated))
})

test_that("E8-05.7: treatment-related attrition stays zero without treatment assignment", {
  panel <- data.frame(
    unit_id = c(
      rep(1L, 4L),
      rep(2L, 4L),
      rep(3L, 3L),
      rep(4L, 3L)
    ),
    year = c(
      1L, 2L, 3L, 4L,
      1L, 2L, 3L, 4L,
      1L, 2L, 3L,
      1L, 2L, 3L
    ),
    y = c(rep(1, 7L), rep(2, 7L)),
    stringsAsFactors = FALSE
  )

  diagnosis <- diagnose_selection_mechanism(
    panel,
    ivar = "unit_id",
    tvar = "year",
    y = "y",
    gvar = NULL,
    verbose = FALSE
  )

  expect_equal(diagnosis$attrition_analysis$overall_attrition, 0.5, tolerance = 1e-12)
  expect_equal(
    diagnosis$attrition_analysis$treatment_related_attrition,
    0,
    tolerance = 1e-12
  )
})

test_that("E8-05.6: all-missing panel returns HIGH risk without crashing", {
  panel <- make_selection_all_missing_panel()

  diagnosis <- diagnose_selection_mechanism(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = "gvar",
    verbose = FALSE
  )

  expect_s3_class(diagnosis, "lwdid_selection_diagnosis")
  expect_equal(diagnosis$missing_rate_overall, 1, tolerance = 1e-12)
  expect_identical(toupper(diagnosis$selection_risk), "HIGH")
  expect_gte(diagnosis$selection_risk_score, 50)
  expect_true(length(diagnosis$warnings) >= 2L)
  expect_true(all(diagnosis$unit_stats$n_observed == 0L))
})

test_that("E8-05 common-timing diagnosis uses d to identify treated units", {
  panel <- make_selection_common_timing_panel()

  diagnosis <- diagnose_selection_mechanism(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = NULL,
    d = "d",
    verbose = FALSE
  )

  treated <- diagnosis$unit_stats[diagnosis$unit_stats$unit_id == 1, ]
  control <- diagnosis$unit_stats[diagnosis$unit_stats$unit_id == 2, ]

  expect_equal(diagnosis$balance_statistics$n_treated_units, 1L)
  expect_equal(diagnosis$balance_statistics$pct_usable_demean, 100, tolerance = 1e-12)
  expect_equal(diagnosis$balance_statistics$pct_usable_detrend, 0, tolerance = 1e-12)
  expect_true(treated$is_treated)
  expect_equal(treated$cohort, 3L)
  expect_equal(treated$n_pre_treatment, 1L)
  expect_equal(treated$n_post_treatment, 1L)
  expect_true(treated$can_use_demean)
  expect_false(treated$can_use_detrend)
  expect_false(control$is_treated)
  expect_true(is.na(control$cohort))
})

test_that("E8-05.6: risk scoring honors strict factor boundaries and warnings", {
  cases <- list(
    list(
      case_id = "mnar_pattern_warning",
      pattern = "MNAR",
      attrition = 0.09,
      dropout_before = 0L,
      dropout_after = 0L,
      balance_ratio = 1.0,
      expected_score = 30,
      expected_risk = "medium",
      expected_warning_count = 1L,
      expected_factors = list(
        missing_pattern = 30,
        attrition = 0,
        differential_attrition = 0,
        balance = 0
      ),
      warning_pattern = "selection on unobservables"
    ),
    list(
      case_id = "attrition_below_10pct",
      pattern = "MCAR",
      attrition = 0.09,
      dropout_before = 0L,
      dropout_after = 0L,
      balance_ratio = 1.0,
      expected_score = 0,
      expected_risk = "low",
      expected_warning_count = 0L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 0,
        differential_attrition = 0,
        balance = 0
      ),
      warning_pattern = NULL
    ),
    list(
      case_id = "attrition_at_10pct",
      pattern = "MCAR",
      attrition = 0.10,
      dropout_before = 0L,
      dropout_after = 0L,
      balance_ratio = 1.0,
      expected_score = 12,
      expected_risk = "low",
      expected_warning_count = 0L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 12,
        differential_attrition = 0,
        balance = 0
      ),
      warning_pattern = NULL
    ),
    list(
      case_id = "attrition_at_30pct",
      pattern = "MCAR",
      attrition = 0.30,
      dropout_before = 0L,
      dropout_after = 0L,
      balance_ratio = 1.0,
      expected_score = 25,
      expected_risk = "medium",
      expected_warning_count = 1L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 25,
        differential_attrition = 0,
        balance = 0
      ),
      warning_pattern = "High attrition rate"
    ),
    list(
      case_id = "differential_attrition_zero_baseline_guard",
      pattern = "MCAR",
      attrition = 0.10,
      dropout_before = 0L,
      dropout_after = 5L,
      balance_ratio = 1.0,
      expected_score = 12,
      expected_risk = "low",
      expected_warning_count = 0L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 12,
        differential_attrition = 0,
        balance = 0
      ),
      warning_pattern = NULL
    ),
    list(
      case_id = "differential_attrition_gt_1_5x",
      pattern = "MCAR",
      attrition = 0.10,
      dropout_before = 2L,
      dropout_after = 4L,
      balance_ratio = 1.0,
      expected_score = 27,
      expected_risk = "medium",
      expected_warning_count = 0L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 12,
        differential_attrition = 15,
        balance = 0
      ),
      warning_pattern = NULL
    ),
    list(
      case_id = "differential_attrition_gt_2x",
      pattern = "MCAR",
      attrition = 0.10,
      dropout_before = 2L,
      dropout_after = 5L,
      balance_ratio = 1.0,
      expected_score = 37,
      expected_risk = "medium",
      expected_warning_count = 1L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 12,
        differential_attrition = 25,
        balance = 0
      ),
      warning_pattern = "dropout after treatment"
    ),
    list(
      case_id = "balance_ratio_at_0_8",
      pattern = "MCAR",
      attrition = 0.10,
      dropout_before = 0L,
      dropout_after = 0L,
      balance_ratio = 0.80,
      expected_score = 22,
      expected_risk = "low",
      expected_warning_count = 0L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 12,
        differential_attrition = 0,
        balance = 10
      ),
      warning_pattern = NULL
    ),
    list(
      case_id = "balance_ratio_above_0_5",
      pattern = "MCAR",
      attrition = 0.10,
      dropout_before = 0L,
      dropout_after = 0L,
      balance_ratio = 0.60,
      expected_score = 22,
      expected_risk = "low",
      expected_warning_count = 0L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 12,
        differential_attrition = 0,
        balance = 10
      ),
      warning_pattern = NULL
    ),
    list(
      case_id = "balance_ratio_at_0_5",
      pattern = "MCAR",
      attrition = 0.10,
      dropout_before = 0L,
      dropout_after = 0L,
      balance_ratio = 0.50,
      expected_score = 32,
      expected_risk = "medium",
      expected_warning_count = 1L,
      expected_factors = list(
        missing_pattern = 0,
        attrition = 12,
        differential_attrition = 0,
        balance = 20
      ),
      warning_pattern = "Low balance ratio"
    )
  )

  for (case in cases) {
    result <- .assess_selection_risk(
      missing_rates = list(
        overall_missing_rate = case$attrition,
        is_balanced = FALSE
      ),
      missing_pattern = list(pattern = case$pattern),
      attrition = list(
        overall_attrition = case$attrition,
        dropout_before_treatment = case$dropout_before,
        dropout_after_treatment = case$dropout_after
      ),
      balance = list(
        is_balanced = FALSE,
        balance_ratio = case$balance_ratio
      )
    )

    expect_equal(
      result$score,
      case$expected_score,
      info = case$case_id
    )
    expect_identical(
      result$risk,
      case$expected_risk,
      info = case$case_id
    )
    expect_equal(
      result$factors,
      case$expected_factors,
      info = case$case_id
    )
    expect_equal(
      length(result$warnings),
      case$expected_warning_count,
      info = case$case_id
    )

    if (!is.null(case$warning_pattern)) {
      expect_true(
        any(grepl(case$warning_pattern, result$warnings, ignore.case = TRUE)),
        info = case$case_id
      )
    }
  }
})

test_that("E8-05.6: single-unit diagnosis keeps balance stats finite and usable", {
  panel <- make_selection_single_unit_panel()

  balance <- .compute_balance_statistics(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = "gvar",
    never_treated_values = c(0, Inf)
  )
  diagnosis <- diagnose_selection_mechanism(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = "gvar",
    verbose = FALSE
  )

  expect_equal(balance$n_units, 1L)
  expect_equal(balance$n_treated_units, 1L)
  expect_equal(balance$std_obs_per_unit, 0, tolerance = 1e-12)
  expect_equal(balance$balance_ratio, 1, tolerance = 1e-12)
  expect_equal(balance$pct_usable_demean, 100, tolerance = 1e-12)
  expect_equal(balance$pct_usable_detrend, 0, tolerance = 1e-12)

  expect_s3_class(diagnosis, "lwdid_selection_diagnosis")
  expect_equal(diagnosis$missing_pattern, "MCAR")
  expect_equal(diagnosis$missing_pattern_confidence, 1, tolerance = 1e-12)
  expect_equal(diagnosis$selection_risk_score, 0)
  expect_identical(toupper(diagnosis$selection_risk), "LOW")
  expect_length(diagnosis$warnings, 0L)
  expect_length(diagnosis$recommendations, 2L)
  expect_true(any(grepl("detrend", diagnosis$recommendations, ignore.case = TRUE)))
  expect_equal(
    diagnosis$balance_statistics$std_obs_per_unit,
    0,
    tolerance = 1e-12
  )
  expect_equal(
    diagnosis$balance_statistics$balance_ratio,
    1,
    tolerance = 1e-12
  )
  expect_equal(nrow(diagnosis$unit_stats), 1L)
  expect_true(diagnosis$unit_stats$is_treated[[1L]])
})

test_that("E8-05.7: covariates end-to-end path runs MAR regression without formula parse errors", {
  panel <- make_selection_covariate_mar_panel()
  diagnosis <- NULL

  expect_no_error({
    diagnosis <- diagnose_selection_mechanism(
      panel,
      ivar = "id",
      tvar = "time",
      y = "y",
      gvar = "gvar",
      covariates = "x",
      verbose = FALSE
    )
  })

  expect_s3_class(diagnosis, "lwdid_selection_diagnosis")
  expect_true(any(vapply(
    diagnosis$selection_tests,
    function(test) identical(test$test_name, "Selection on Observables (MAR) Test"),
    logical(1)
  )))
  expect_true(diagnosis$selection_risk %in% c("low", "medium", "high"))
  expect_true(is.finite(diagnosis$selection_risk_score))
  expect_gte(diagnosis$selection_risk_score, 0)
  expect_lte(diagnosis$selection_risk_score, 100)
  expect_equal(diagnosis$missing_rate_overall, 4 / 30, tolerance = 1e-12)
  expect_equal(diagnosis$attrition_analysis$n_partial, 3L)
})

test_that("E8-05.7: non-finite MCAR/MNAR p-values do not crash diagnosis", {
  panel <- do.call(
    rbind,
    lapply(seq_len(12L), function(unit_id) {
      years <- if (unit_id <= 6L) {
        1:4
      } else {
        1:3
      }
      data.frame(
        id = unit_id,
        time = years,
        y = 1,
        stringsAsFactors = FALSE
      )
    })
  )

  diagnosis <- NULL
  expect_no_error({
    diagnosis <- diagnose_selection_mechanism(
      panel,
      ivar = "id",
      tvar = "time",
      y = "y",
      verbose = FALSE
    )
  })

  expect_s3_class(diagnosis, "lwdid_selection_diagnosis")
  expect_true(is.finite(diagnosis$selection_risk_score))
  expect_false(identical(toupper(diagnosis$missing_pattern), "MNAR"))

  if (length(diagnosis$selection_tests) > 0L) {
    reject_flags <- vapply(
      diagnosis$selection_tests,
      function(test) isTRUE(test$reject_null),
      logical(1)
    )
    expect_false(any(reject_flags))
  }
})

test_that("E8-05.7: layer5 attrition surface exposes cohort ordering and accounting identity", {
  panel <- make_selection_layer5_attrition_panel()

  diagnosis <- diagnose_selection_mechanism(
    panel,
    ivar = "id",
    tvar = "time",
    y = "y",
    gvar = "gvar",
    verbose = FALSE
  )

  expect_s3_class(diagnosis, "lwdid_selection_diagnosis")
  expect_type(diagnosis$attrition_analysis$attrition_by_cohort, "list")
  expect_equal(
    diagnosis$attrition_analysis$attrition_by_cohort[["6"]],
    0.5,
    tolerance = 1e-12
  )
  expect_equal(
    diagnosis$attrition_analysis$attrition_by_cohort[["8"]],
    0.25,
    tolerance = 1e-12
  )
  expect_gt(
    diagnosis$attrition_analysis$attrition_by_cohort[["6"]],
    diagnosis$attrition_analysis$attrition_by_cohort[["8"]]
  )

  expect_type(diagnosis$attrition_analysis$attrition_by_period, "list")
  expect_equal(
    diagnosis$attrition_analysis$attrition_by_period[["10"]],
    4 / 12,
    tolerance = 1e-12
  )
  expect_equal(
    diagnosis$attrition_by_period[["10"]],
    4 / 12,
    tolerance = 1e-12
  )
  expect_gt(
    diagnosis$attrition_analysis$attrition_by_period[["10"]],
    diagnosis$attrition_analysis$attrition_by_period[["9"]]
  )

  expect_equal(diagnosis$attrition_analysis$dropout_before_treatment, 0L)
  expect_equal(diagnosis$attrition_analysis$dropout_after_treatment, 3L)
  expect_equal(sum(diagnosis$unit_stats$is_treated), 8L)
  expect_equal(
    sort(unique(diagnosis$unit_stats$cohort[diagnosis$unit_stats$is_treated])),
    c(6L, 8L)
  )

  expect_equal(
    sum(diagnosis$unit_stats$n_observed + diagnosis$unit_stats$n_missing),
    diagnosis$balance_statistics$n_units * diagnosis$balance_statistics$n_periods
  )
  expect_equal(nrow(diagnosis$unit_stats), diagnosis$balance_statistics$n_units)
})

test_that("E8-05.7: unit_id-named ivar preserves mixed cohort treatment flags", {
  panel <- make_selection_layer5_attrition_panel()
  names(panel)[names(panel) == "id"] <- "unit_id"
  names(panel)[names(panel) == "time"] <- "year"

  diagnosis <- diagnose_selection_mechanism(
    panel,
    ivar = "unit_id",
    tvar = "year",
    y = "y",
    gvar = "gvar",
    verbose = FALSE
  )

  expect_equal(sum(diagnosis$unit_stats$is_treated), 8L)
  expect_equal(
    sort(unique(diagnosis$unit_stats$cohort[diagnosis$unit_stats$is_treated])),
    c(6L, 8L)
  )
})
