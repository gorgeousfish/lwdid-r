# ============================================================================
# test-advanced-usage-doc-contract.R
# Advanced vignette estimator examples must match real Castle execution.
# ============================================================================

test_that("advanced estimator docs spell out the cohort control group", {
  doc_paths <- c(
    resolve_package_source_file("vignettes", "advanced-usage.Rmd"),
    resolve_package_source_file("inst", "doc", "advanced-usage.Rmd"),
    resolve_package_source_file("doc", "advanced-usage.Rmd")
  )

  for (path in doc_paths) {
    expect_true(file.exists(path), info = sprintf("Missing doc file: %s", path))
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")

    for (estimator in c("ipw", "ipwra", "psm")) {
      pattern <- paste0(
        'estimator = "', estimator, '"[\\s\\S]*',
        'aggregate = "cohort", control_group = "never_treated"'
      )
      expect_match(
        text,
        pattern,
        perl = TRUE,
        info = sprintf("%s must pin %s cohort examples to never-treated controls", path, estimator)
      )
    }
  }
})

test_that("advanced estimator docs state propensity-score controls explicitly", {
  doc_paths <- c(
    resolve_package_source_file("README.md"),
    resolve_package_source_file("vignettes", "advanced-usage.Rmd"),
    resolve_package_source_file("inst", "doc", "advanced-usage.Rmd"),
    resolve_package_source_file("doc", "advanced-usage.Rmd")
  )

  for (path in doc_paths) {
    expect_true(file.exists(path), info = sprintf("Missing doc file: %s", path))
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")

    expect_match(
      text,
      "ps_controls",
      fixed = TRUE,
      info = sprintf("%s should name explicit propensity-score controls", path)
    )
    expect_match(
      text,
      "fall back\\s+to\\s+`controls`",
      perl = TRUE,
      info = sprintf("%s should document propensity-control fallback semantics", path)
    )
  }

  readme_text <- paste(
    readLines(resolve_package_source_file("README.md"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    readme_text,
    'estimator = "ipw"[\\s\\S]*ps_controls = c\\("police", "unemployrt", "income"\\)',
    perl = TRUE,
    info = "README IPW example should use explicit propensity-score controls"
  )

  vignette_text <- paste(
    readLines(resolve_package_source_file("vignettes", "advanced-usage.Rmd"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    vignette_text,
    'estimator = "ipw"[\\s\\S]*ps_controls = c\\("police", "income", "poverty"\\)',
    perl = TRUE
  )
  expect_match(
    vignette_text,
    'estimator = "psm"[\\s\\S]*ps_controls = c\\("police", "income", "poverty"\\)',
    perl = TRUE
  )
})

test_that("advanced estimator examples run on Castle without cohort-control rewrites", {
  data("castle", package = "lwdid")

  run_example <- function(estimator) {
    warnings <- character()
    result <- withCallingHandlers(
      lwdid(
        data = castle,
        y = "lhomicide",
        ivar = "sid",
        tvar = "year",
        gvar = "gvar",
        rolling = "demean",
        estimator = estimator,
        controls = if (estimator == "ipwra") c("police", "income", "poverty") else NULL,
        ps_controls = c("police", "income", "poverty"),
        aggregate = "cohort",
        control_group = "never_treated"
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(result = result, warnings = warnings)
  }

  for (estimator in c("ipw", "ipwra", "psm")) {
    example <- run_example(estimator)

    expect_s3_class(example$result, "lwdid_result")
    expect_true(is.finite(example$result$att))
    expect_false(
      any(grepl("control_group.*never_treated|not_yet_treated", example$warnings)),
      info = sprintf("%s example should not rely on automatic control-group rewriting", estimator)
    )
  }
})

test_that("advanced parallel and RI examples match real-data execution", {
  data("castle", package = "lwdid")
  data("smoking", package = "lwdid")

  doc_text <- paste(
    readLines(
      resolve_package_source_file("vignettes", "advanced-usage.Rmd"),
      warn = FALSE
    ),
    collapse = "\n"
  )

  expect_match(
    doc_text,
    'ri = TRUE, ri_method = "permutation"',
    info = "RI example should use the stable permutation method on smoking data"
  )
  expect_match(
    doc_text,
    'control_group = "never_treated", parallel = TRUE',
    info = "parallel cohort example should pin the cross-cohort control group"
  )
  expect_false(
    grepl("Exact p-values", doc_text, fixed = TRUE),
    info = "advanced RI prose should not overstate simulated randomization inference"
  )

  ri_warnings <- character()
  ri_result <- withCallingHandlers(
    lwdid(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      ri = TRUE,
      ri_method = "permutation",
      rireps = 100L,
      seed = 42
    ),
    warning = function(w) {
      ri_warnings <<- c(ri_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(ri_result, "lwdid_result")
  expect_true(is.finite(ri_result$ri_pvalue))
  expect_equal(ri_result$ri_method, "permutation")
  expect_equal(ri_result$ri_valid, 100L)
  expect_false(any(grepl("replications failed|Randomization inference failed", ri_warnings)))

  parallel_warnings <- character()
  parallel_result <- withCallingHandlers(
    lwdid(
      data = castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      aggregate = "cohort",
      control_group = "never_treated",
      parallel = TRUE,
      n_cores = 2
    ),
    warning = function(w) {
      parallel_warnings <<- c(parallel_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(parallel_result, "lwdid_result")
  expect_true(is.finite(parallel_result$att))
  expect_false(any(grepl("not_yet_treated|auto-switching", parallel_warnings)))
})

test_that("RI docs describe staggered target statistic semantics", {
  doc_paths <- c(
    resolve_package_source_file("README.md"),
    resolve_package_source_file("vignettes", "advanced-usage.Rmd"),
    resolve_package_source_file("inst", "doc", "advanced-usage.Rmd"),
    resolve_package_source_file("doc", "advanced-usage.Rmd")
  )

  for (path in doc_paths) {
    expect_true(file.exists(path), info = sprintf("Missing doc file: %s", path))
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")

    expect_match(
      text,
      "ri_observed_stat",
      fixed = TRUE,
      info = sprintf("%s should identify the observed RI statistic field", path)
    )
    expect_match(
      text,
      "ri_target",
      fixed = TRUE,
      info = sprintf("%s should identify the RI target metadata field", path)
    )
    expect_match(
      text,
      'aggregate = "none"[\\s\\S]*first finite cohort-time ATT',
      perl = TRUE,
      info = sprintf("%s should state the staggered aggregate='none' RI target", path)
    )
  }

  help_text <- paste(
    readLines(resolve_package_source_file("R", "lwdid.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    help_text,
    'ri_pvalue[\\s\\S]*ri_observed_stat[\\s\\S]*aggregate = "none"[\\s\\S]*first finite cohort-time ATT',
    perl = TRUE,
    info = "lwdid() roxygen should document the RI observed-statistic boundary"
  )
})

test_that("WCB examples request reproducible native bootstrap execution", {
  data("smoking", package = "lwdid")

  doc_paths <- c(
    resolve_package_source_file("README.md"),
    resolve_package_source_file("vignettes", "advanced-usage.Rmd"),
    resolve_package_source_file("inst", "doc", "advanced-usage.Rmd"),
    resolve_package_source_file("doc", "advanced-usage.Rmd")
  )

  for (path in doc_paths) {
    expect_true(file.exists(path), info = sprintf("Missing doc file: %s", path))
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")

    expect_match(
      text,
      'vce = "bootstrap"[\\s\\S]*wcb_seed = 42[\\s\\S]*use_fwildclusterboot = FALSE',
      perl = TRUE,
      info = sprintf("%s should document reproducible native WCB execution", path)
    )
    expect_match(
      text,
      "requested and actual bootstrap draws",
      fixed = TRUE,
      info = sprintf("%s should explain the WCB requested/actual draw contract", path)
    )
    expect_match(
      text,
      "full sign-pattern\\s+enumeration",
      perl = TRUE,
      info = sprintf("%s should explain why WCB actual draws can differ", path)
    )
  }

  readme_text <- paste(
    readLines(resolve_package_source_file("README.md"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(
    readme_text,
    "automatically upgrades to WCB only for fewer than 20 clusters",
    fixed = TRUE
  )
  expect_false(
    grepl('vce = "cluster"[\\s\\S]*aggregate = "overall"', readme_text),
    info = "README WCB example should not imply Castle cluster VCE is WCB"
  )

  wcb_warnings <- character()
  wcb_messages <- character()
  wcb_result <- withCallingHandlers(
    lwdid(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      rolling = "demean",
      vce = "bootstrap",
      cluster_var = "state",
      wcb_reps = 99L,
      wcb_seed = 42L,
      use_fwildclusterboot = FALSE
    ),
    warning = function(w) {
      wcb_warnings <<- c(wcb_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      wcb_messages <<- c(wcb_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_s3_class(wcb_result, "lwdid_result")
  expect_equal(wcb_result$vce_type, "wild_cluster_bootstrap")
  expect_true(is.finite(wcb_result$pvalue_wcb))
  expect_equal(wcb_result$wcb_details$requested_n_bootstrap, 99L)
  expect_equal(wcb_result$wcb_details$actual_n_bootstrap, 99L)
  expect_equal(wcb_result$wcb_details$n_clusters, 39L)
  expect_false(any(grepl("fwildclusterboot not installed", wcb_messages)))
  expect_false(any(grepl("auto-enabling Wild Cluster Bootstrap", wcb_messages)))
  expect_false(any(grepl("replications failed", wcb_warnings)))
})
