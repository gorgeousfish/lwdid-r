# ============================================================================
# test-staggered-overall-export-contract.R
# Overall staggered exports must preserve aggregation weights from real data.
# ============================================================================

test_that("overall aggregation docs explain exported weight semantics", {
  doc_paths <- c(
    resolve_package_source_file("README.md"),
    resolve_package_source_file("vignettes", "getting-started.Rmd"),
    resolve_package_source_file("inst", "doc", "getting-started.Rmd")
  )

  for (path in doc_paths) {
    expect_true(file.exists(path), info = sprintf("Missing doc file: %s", path))
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")

    expect_match(
      text,
      "treated-unit cohort sizes",
      fixed = TRUE,
      info = sprintf("%s must explain the overall aggregation weights", path)
    )
    expect_match(
      text,
      "cohort_weights",
      fixed = TRUE,
      info = sprintf("%s must document exported cohort weights", path)
    )
    expect_match(
      text,
      "effective_weights",
      fixed = TRUE,
      info = sprintf("%s must document effective weight diagnostics", path)
    )
  }
})

test_that("Castle overall exports preserve cohort aggregation weights", {
  data("castle", package = "lwdid")

  result <- suppressWarnings(lwdid(
    data = castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    aggregate = "overall",
    control_group = "never_treated"
  ))

  expect_s3_class(result, "lwdid_result")
  expect_equal(result$att, result$overall_effect$att, tolerance = 1e-12)
  expect_equal(result$n_units, length(unique(castle$sid)))
  expect_equal(result$n_periods, length(unique(castle$year)))
  expect_equal(result$n_cohorts, length(unique(result$att_by_cohort$cohort)))
  expect_false(is.null(result$cohort_weights))
  expect_false(is.null(result$effective_weights))
  expect_equal(
    unlist(result$cohort_weights),
    unlist(result$overall_effect$cohort_weights),
    tolerance = 1e-12
  )
  expect_equal(
    unlist(result$effective_weights),
    unlist(result$overall_effect$effective_weights),
    tolerance = 1e-12
  )
  expect_equal(sum(unlist(result$cohort_weights)), 1, tolerance = 1e-12)

  expected_castle_weights <- c(
    "2005" = 1 / 21,
    "2006" = 13 / 21,
    "2007" = 4 / 21,
    "2008" = 2 / 21,
    "2009" = 1 / 21
  )
  expect_equal(
    unlist(result$cohort_weights)[names(expected_castle_weights)],
    expected_castle_weights,
    tolerance = 1e-12
  )

  by_cohort_file <- tempfile(fileext = ".csv")
  on.exit(unlink(by_cohort_file), add = TRUE)
  by_cohort <- suppressMessages(to_csv(
    result,
    file = by_cohort_file,
    what = "by_cohort"
  ))

  expected_weights <- unlist(result$cohort_weights)
  expected_sizes <- result$cohort_sample_sizes[as.character(by_cohort$cohort)]

  expect_true("weight" %in% names(by_cohort))
  expect_true("effective_weight" %in% names(by_cohort))
  expect_true("n" %in% names(by_cohort))
  expect_equal(
    by_cohort$weight,
    as.numeric(expected_weights[as.character(by_cohort$cohort)]),
    tolerance = 1e-12
  )
  expect_equal(by_cohort$n, as.integer(expected_sizes))
  expect_equal(
    by_cohort$effective_weight,
    as.numeric(unlist(result$effective_weights)[as.character(by_cohort$cohort)]),
    tolerance = 1e-12
  )

  summary_file <- tempfile(fileext = ".csv")
  on.exit(unlink(summary_file), add = TRUE)
  summary_df <- suppressMessages(to_csv(
    result,
    file = summary_file,
    what = "summary"
  ))
  expect_equal(summary_df$n_cohorts, length(unique(result$att_by_cohort$cohort)))

  dict <- to_dict(result)
  expect_equal(dict$n_cohorts, length(unique(result$att_by_cohort$cohort)))
  expect_equal(
    unlist(dict$cohort_weights),
    unlist(result$cohort_weights),
    tolerance = 1e-12
  )
  expect_equal(
    unlist(dict$effective_weights),
    unlist(result$effective_weights),
    tolerance = 1e-12
  )
})

test_that("Castle non-RA overall aggregates non-RA cohort effects", {
  data("castle", package = "lwdid")

  result <- suppressWarnings(lwdid(
    data = castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ipw",
    aggregate = "overall",
    control_group = "never_treated",
    controls = c("income", "unemployrt", "poverty")
  ))

  expect_s3_class(result, "lwdid_result")
  expect_equal(result$estimator, "ipw")
  expect_false(is.null(result$overall_effect))
  expect_equal(
    result$overall_effect$aggregation_method,
    "cohort_size_weighted_non_ra_cohort_effects"
  )
  expect_equal(
    result$overall_effect$se_aggregation,
    "diagonal_weighted_cohort_se"
  )
  expect_equal(
    result$overall_effect$covariance_assumption,
    "zero_cross_cohort_covariance"
  )

  cohort_df <- result$att_by_cohort[
    is.finite(result$att_by_cohort$att) & is.finite(result$att_by_cohort$se), ,
    drop = FALSE
  ]
  weights <- unlist(result$cohort_weights)
  cohort_weights <- as.numeric(weights[as.character(cohort_df$cohort)])
  expected_att <- sum(cohort_weights * cohort_df$att)
  expected_se <- sqrt(sum(cohort_weights^2 * cohort_df$se^2))
  expected_p <- 2 * stats::pnorm(abs(expected_att / expected_se), lower.tail = FALSE)

  expect_equal(result$att_overall, expected_att, tolerance = 1e-12)
  expect_equal(result$se_overall, expected_se, tolerance = 1e-12)
  expect_equal(result$pvalue_overall, expected_p, tolerance = 1e-12)
  expect_equal(result$att, result$att_overall, tolerance = 1e-12)
  expect_equal(result$se_att, result$se_overall, tolerance = 1e-12)
})

test_that("Castle overall exports distinguish effective weights after regression-row drops", {
  data("castle", package = "lwdid")

  castle_missing <- castle
  drop_ids <- unique(castle_missing$sid[
    !is.na(castle_missing$gvar) & castle_missing$gvar == 2006
  ])[1:3]
  castle_missing$income[castle_missing$sid %in% drop_ids] <- NA_real_

  warnings <- character()
  result <- withCallingHandlers(
    lwdid(
      data = castle_missing,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      estimator = "ra",
      aggregate = "overall",
      control_group = "never_treated",
      controls = "income"
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expected_computation_weights <- c(
    "2005" = 1 / 21,
    "2006" = 13 / 21,
    "2007" = 4 / 21,
    "2008" = 2 / 21,
    "2009" = 1 / 21
  )
  expected_effective_weights <- c(
    "2005" = 1 / 18,
    "2006" = 10 / 18,
    "2007" = 4 / 18,
    "2008" = 2 / 18,
    "2009" = 1 / 18
  )

  expect_equal(
    unlist(result$cohort_weights)[names(expected_computation_weights)],
    expected_computation_weights,
    tolerance = 1e-12
  )
  expect_equal(
    unlist(result$effective_weights)[names(expected_effective_weights)],
    expected_effective_weights,
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    unlist(result$cohort_weights),
    unlist(result$effective_weights),
    tolerance = 1e-12
  )))
  expect_true(any(grepl("Post-dropna cohort weight deviation", warnings)))

  by_cohort_file <- tempfile(fileext = ".csv")
  on.exit(unlink(by_cohort_file), add = TRUE)
  by_cohort <- suppressMessages(to_csv(
    result,
    file = by_cohort_file,
    what = "by_cohort"
  ))

  expect_true("weight" %in% names(by_cohort))
  expect_true("effective_weight" %in% names(by_cohort))
  expect_equal(
    by_cohort$weight,
    as.numeric(expected_computation_weights[as.character(by_cohort$cohort)]),
    tolerance = 1e-12
  )
  expect_equal(
    by_cohort$effective_weight,
    as.numeric(expected_effective_weights[as.character(by_cohort$cohort)]),
    tolerance = 1e-12
  )

  dict <- to_dict(result)
  expect_equal(
    unlist(dict$cohort_weights)[names(expected_computation_weights)],
    expected_computation_weights,
    tolerance = 1e-12
  )
  expect_equal(
    unlist(dict$effective_weights)[names(expected_effective_weights)],
    expected_effective_weights,
    tolerance = 1e-12
  )
})
