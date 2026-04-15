test_that("E8-07.1: get_diagnostics returns normalized diagnostics entries", {
  selection_diag <- structure(
    list(selection_risk = "Low", attrition_rate = 0.05),
    class = "lwdid_selection_diagnosis"
  )
  trend_diag <- structure(
    list(recommended_method = "demean", confidence_level = "High"),
    class = "lwdid_transformation_recommendation"
  )

  result <- new_lwdid_result(
    att = 1,
    se_att = 0.2,
    df_inference = 10L,
    diagnostics = list(
      selection = selection_diag,
      parallel_trends = trend_diag
    )
  )

  all_diags <- get_diagnostics(result, "all")

  expect_named(
    all_diags,
    c("clustering", "selection", "trends", "sensitivity")
  )
  expect_identical(all_diags$selection, selection_diag)
  expect_identical(all_diags$trends, trend_diag)
  expect_null(all_diags$clustering)
  expect_null(all_diags$sensitivity)
})

test_that("E8-07.1: get_diagnostics returns type-specific diagnostics", {
  selection_diag <- structure(
    list(selection_risk = "Low", attrition_rate = 0.05),
    class = "lwdid_selection_diagnosis"
  )

  result <- new_lwdid_result(
    att = 1,
    se_att = 0.2,
    df_inference = 10L,
    diagnostics = list(selection = selection_diag)
  )

  expect_identical(get_diagnostics(result, "selection"), selection_diag)
  expect_null(get_diagnostics(result, "sensitivity"))
})

test_that("E8-07.1: get_diagnostics messages when diagnostics were not run", {
  result <- new_lwdid_result(
    att = 1,
    se_att = 0.2,
    df_inference = 10L,
    diagnostics = list(controls_tier = "none")
  )

  expect_message(
    returned <- get_diagnostics(result),
    "return_diagnostics=TRUE"
  )
  expect_null(returned)
})

test_that("E8-07.1: get_diagnostics validates input object and type", {
  result <- new_lwdid_result(att = 1, se_att = 0.2, df_inference = 10L)

  expect_error(get_diagnostics(list(), "all"), "lwdid")
  expect_error(get_diagnostics(result, "invalid"), "all, clustering, selection, trends, sensitivity")
})

test_that("E8-07.2: lwdid_diagnose matches direct diagnostics on castle real data", {
  data(castle, package = "lwdid")

  direct_clustering <- suppressWarnings(diagnose_clustering(
    castle,
    ivar = "sid",
    potential_cluster_vars = c("state_name"),
    gvar = "gvar",
    verbose = FALSE
  ))
  direct_selection <- suppressWarnings(diagnose_selection_mechanism(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    covariates = c("income", "unemployrt", "poverty"),
    verbose = FALSE
  ))
  direct_trends <- suppressWarnings(lwdid_recommend_transformation(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    run_all_diagnostics = FALSE,
    verbose = FALSE
  ))

  suite <- suppressWarnings(lwdid_diagnose(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_vars = c("state_name"),
    verbose = FALSE
  ))

  expect_s3_class(suite, "lwdid_diagnostics_suite")
  expect_identical(
    suite$clustering$recommended_cluster_var,
    direct_clustering$recommended_cluster_var
  )
  expect_identical(
    suite$clustering$treatment_variation_level,
    direct_clustering$treatment_variation_level
  )
  expect_identical(
    suite$selection$selection_risk,
    direct_selection$selection_risk
  )
  expect_equal(
    suite$selection$missing_rate_overall,
    direct_selection$missing_rate_overall
  )
  expect_identical(
    suite$trends$recommended_method,
    direct_trends$recommended_method
  )
  expect_identical(
    suite$trends$confidence_level,
    direct_trends$confidence_level
  )
})

test_that("E8-07.2: lwdid_diagnose handles common-timing smoking data", {
  data(smoking, package = "lwdid")

  direct_selection <- diagnose_selection_mechanism(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    covariates = c("lretprice", "lnincome"),
    verbose = FALSE
  )
  direct_trends <- lwdid_recommend_transformation(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    gvar = NULL,
    controls = c("lretprice", "lnincome"),
    run_all_diagnostics = FALSE,
    verbose = FALSE
  )

  suite <- lwdid_diagnose(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    controls = c("lretprice", "lnincome"),
    cluster_vars = NULL,
    verbose = FALSE
  )

  expect_s3_class(suite, "lwdid_diagnostics_suite")
  expect_null(suite$clustering)
  expect_identical(
    suite$selection$selection_risk,
    direct_selection$selection_risk
  )
  expect_equal(
    suite$selection$missing_rate_overall,
    direct_selection$missing_rate_overall
  )
  expect_identical(
    suite$trends$recommended_method,
    direct_trends$recommended_method
  )
  expect_identical(
    suite$trends$confidence_level,
    direct_trends$confidence_level
  )
})

test_that("E8-07.2: lwdid_diagnose respects custom never-treated markers", {
  data(castle, package = "lwdid")

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

  expect_s3_class(custom_suite, "lwdid_diagnostics_suite")
  expect_null(custom_suite$clustering)
  expect_identical(
    custom_suite$selection$selection_risk,
    zero_suite$selection$selection_risk
  )
  expect_equal(
    custom_suite$selection$missing_rate_overall,
    zero_suite$selection$missing_rate_overall
  )
  expect_identical(
    custom_suite$trends$recommended_method,
    zero_suite$trends$recommended_method
  )
  expect_identical(
    custom_suite$trends$confidence_level,
    zero_suite$trends$confidence_level
  )
})

test_that("E8-07.2: lwdid_diagnose isolates failing modules", {
  data(castle, package = "lwdid")

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

  expect_s3_class(suite, "lwdid_diagnostics_suite")
  expect_null(suite$clustering)
  expect_s3_class(suite$selection, "lwdid_selection_diagnosis")
  expect_s3_class(suite$trends, "lwdid_transformation_recommendation")
  expect_true(any(grepl("missing_cluster", caught_warnings, fixed = TRUE)))
})

test_that("E8-07.2: diagnostics suite print and summary expose module summaries", {
  data(castle, package = "lwdid")

  suite <- suppressWarnings(lwdid_diagnose(
    castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    cluster_vars = c("state_name"),
    verbose = FALSE
  ))

  expect_output(print(suite), "^lwdid", perl = TRUE)
  expect_output(summary(suite), "Transformation Recommendation")
})

test_that("E8-07.3: return_diagnostics attaches staggered diagnostics without changing estimates", {
  data(castle, package = "lwdid")

  base_result <- suppressWarnings(lwdid(
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

  diagnosed_result <- suppressWarnings(lwdid(
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

  expect_identical(diagnosed_result$att, base_result$att)
  expect_identical(diagnosed_result$se_att, base_result$se_att)
  expect_s3_class(
    get_diagnostics(diagnosed_result, "selection"),
    "lwdid_selection_diagnosis"
  )
  expect_s3_class(
    get_diagnostics(diagnosed_result, "trends"),
    "lwdid_transformation_recommendation"
  )
  expect_s3_class(
    get_diagnostics(diagnosed_result, "clustering"),
    "lwdid_clustering_diagnosis"
  )
})

test_that("E8-07.3: common-timing mainline diagnostics skip clustering and stay quiet", {
  data(smoking, package = "lwdid")

  base_result <- suppressWarnings(lwdid(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    post = "post",
    verbose = "quiet"
  ))

  captured_messages <- character(0)
  captured_output <- capture.output(
    diagnosed_result <- withCallingHandlers(
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
        captured_messages <<- c(captured_messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )

  expect_identical(diagnosed_result$att, base_result$att)
  expect_identical(diagnosed_result$se_att, base_result$se_att)
  expect_length(captured_output, 0L)
  expect_length(captured_messages, 0L)
  expect_s3_class(
    get_diagnostics(diagnosed_result, "selection"),
    "lwdid_selection_diagnosis"
  )
  expect_s3_class(
    get_diagnostics(diagnosed_result, "trends"),
    "lwdid_transformation_recommendation"
  )
  expect_null(get_diagnostics(diagnosed_result, "clustering"))
})

test_that("E8-07.3: diagnostics plot consumes mainline selection diagnostics shape", {
  skip_if_not_installed("ggplot2")
  data(smoking, package = "lwdid")

  diagnosed_result <- suppressWarnings(lwdid(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    post = "post",
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  selection_plot <- plot(
    diagnosed_result,
    type = "diagnostics",
    which = "selection"
  )

  expect_s3_class(selection_plot, "ggplot")
})

resolve_diagnostics_package_root <- function() {
  if (exists("resolve_package_source_root", mode = "function")) {
    return(resolve_package_source_root())
  }
  normalizePath(".", winslash = "/", mustWork = FALSE)
}

read_diagnostics_namespace <- function() {
  pkg_root <- resolve_diagnostics_package_root()
  parseNamespaceFile(basename(pkg_root), dirname(pkg_root))
}

has_diagnostics_s3_method <- function(ns_info, generic, class_name) {
  s3_methods <- ns_info$S3methods
  if (is.null(s3_methods) || length(s3_methods) == 0L) {
    return(FALSE)
  }

  any(s3_methods[, 1] == generic & s3_methods[, 2] == class_name)
}

test_that("E8-07 namespace contract exports diagnostics integration public surface", {
  ns_info <- read_diagnostics_namespace()

  expect_true("get_diagnostics" %in% ns_info$exports)
  expect_true("lwdid_diagnose" %in% ns_info$exports)
  expect_true(has_diagnostics_s3_method(
    ns_info,
    "print",
    "lwdid_diagnostics_suite"
  ))
  expect_true(has_diagnostics_s3_method(
    ns_info,
    "summary",
    "lwdid_diagnostics_suite"
  ))
})

test_that("E8-07 docs contract includes diagnostics integration Rd assets", {
  pkg_root <- resolve_diagnostics_package_root()
  expected_rd <- c(
    "man/get_diagnostics.Rd",
    "man/lwdid_diagnose.Rd",
    "man/print.lwdid_diagnostics_suite.Rd",
    "man/summary.lwdid_diagnostics_suite.Rd"
  )

  missing_rd <- expected_rd[
    !file.exists(file.path(pkg_root, expected_rd))
  ]

  expect(
    length(missing_rd) == 0L,
    paste("Missing Rd files:", paste(missing_rd, collapse = ", "))
  )
})

test_that("E8-07.5.4: lwdid_diagnose returns all-NULL suite when every module fails", {
  tiny_df <- data.frame(
    id = rep(1:3, each = 2),
    year = rep(1:2, times = 3),
    y = c(1, 2, 3, 4, 5, 6),
    gvar = c(2, 2, 2, 2, Inf, Inf),
    state_name = rep(c("A", "B", "C"), each = 2)
  )

  testthat::local_mocked_bindings(
    diagnose_selection_mechanism = function(...) stop("selection boom"),
    lwdid_recommend_transformation = function(...) stop("trend boom"),
    diagnose_clustering = function(...) stop("cluster boom"),
    .package = "lwdid"
  )

  warnings <- character(0)
  suite <- withCallingHandlers(
    lwdid_diagnose(
      tiny_df,
      y = "y",
      ivar = "id",
      tvar = "year",
      gvar = "gvar",
      cluster_vars = "state_name",
      verbose = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(suite, "lwdid_diagnostics_suite")
  expect_null(suite$selection)
  expect_null(suite$trends)
  expect_null(suite$clustering)
  expect_true(any(grepl("selection boom", warnings, fixed = TRUE)))
  expect_true(any(grepl("trend boom", warnings, fixed = TRUE)))
  expect_true(any(grepl("cluster boom", warnings, fixed = TRUE)))
  suite_output <- paste(capture.output(print(suite)), collapse = "\n")
  expect_true(grepl("Not run", suite_output, fixed = TRUE))
})

test_that("E8-07.5.4: lwdid_diagnose survives tiny-sample input", {
  tiny_df <- data.frame(
    id = rep(1:3, each = 2),
    year = rep(1:2, times = 3),
    y = c(1, 2, 3, 4, 5, 6),
    gvar = c(2, 2, 2, 2, Inf, Inf),
    x = c(0, 1, 1, 0, 0, 1)
  )

  expect_no_error(
    suite <- suppressWarnings(lwdid_diagnose(
      tiny_df,
      y = "y",
      ivar = "id",
      tvar = "year",
      gvar = "gvar",
      controls = "x",
      cluster_vars = NULL,
      verbose = FALSE
    ))
  )
  expect_s3_class(suite, "lwdid_diagnostics_suite")
})

test_that("E8-07.5.4: default mainline keeps only legacy controls metadata when diagnostics are off", {
  data(smoking, package = "lwdid")

  result <- suppressWarnings(lwdid(
    smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    post = "post",
    verbose = "quiet"
  ))

  expect_named(result$diagnostics, "controls_tier")
  expect_true(is.character(result$diagnostics$controls_tier))
  expect_message(returned <- get_diagnostics(result), "return_diagnostics=TRUE")
  expect_null(returned)
})
