# test-lwdid-wcb-integration.R

local_suppress_lwdid_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (inherits(w, "lwdid_warning") && !inherits(w, "lwdid_data")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

capture_conditions <- function(expr) {
  messages <- character(0)
  warnings <- character(0)

  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, messages = messages, warnings = warnings)
}

make_e9_05_common_panel <- function(
    n_units = 12L,
    n_periods = 8L,
    n_treated = 6L,
    n_clusters = 6L,
    seed = 42L
) {
  set.seed(seed)

  ids <- seq_len(n_units)
  years <- 2000L + seq_len(n_periods) - 1L

  panel <- expand.grid(
    id = ids,
    year = years,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  panel <- panel[order(panel$id, panel$year), ]

  unit_x1 <- stats::rnorm(n_units)
  unit_cluster <- rep(seq_len(n_clusters), length.out = n_units)
  unit_effect <- stats::rnorm(n_units, sd = 0.3)

  panel$d <- as.integer(panel$id <= n_treated)
  panel$post <- as.integer(panel$year >= 2004L)
  panel$x1 <- unit_x1[panel$id]
  panel$state <- sprintf("g%02d", unit_cluster[panel$id])

  time_index <- panel$year - min(panel$year)
  panel$y <- 1.5 +
    0.15 * time_index +
    1.25 * panel$d * panel$post +
    0.4 * panel$x1 +
    unit_effect[panel$id] +
    stats::rnorm(nrow(panel), sd = 0.1)

  rownames(panel) <- NULL
  panel
}

make_e9_05_seasonal_panel <- function(
    n_units = 8L,
    n_periods = 12L,
    n_treated = 4L,
    seed = 123L
) {
  set.seed(seed)

  ids <- seq_len(n_units)
  times <- seq_len(n_periods)

  panel <- expand.grid(
    id = ids,
    time = times,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  panel <- panel[order(panel$id, panel$time), ]

  quarter <- ((panel$time - 1L) %% 4L) + 1L
  unit_effect <- stats::rnorm(n_units, sd = 0.2)
  season_shift <- c(0, 0.8, -0.4, 0.6)

  panel$d <- as.integer(panel$id <= n_treated)
  panel$post <- as.integer(panel$time >= 9L)
  panel$quarter <- quarter
  panel$y <- 2 +
    season_shift[quarter] +
    0.5 * panel$d * panel$post +
    unit_effect[panel$id] +
    stats::rnorm(nrow(panel), sd = 0.05)

  rownames(panel) <- NULL
  panel
}

test_that("TC-9.5.1: vce='bootstrap' triggers WCB and attaches integration fields", {
  panel <- make_e9_05_common_panel()
  result <- NULL

  expect_no_error(
    result <- local_suppress_lwdid_warnings(lwdid(
      panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = "x1",
      vce = "bootstrap",
      cluster_var = "state",
      wcb_reps = 39L,
      wcb_type = "mammen",
      wcb_seed = 123L,
      wcb_restricted_model = "intercept_only",
      use_fwildclusterboot = FALSE
    ))
  )

  expect_identical(result$vce_type, "wild_cluster_bootstrap")
  expect_false(is.null(result$pvalue_wcb))
  expect_false(is.null(result$ci_lower_wcb))
  expect_false(is.null(result$ci_upper_wcb))
  expect_false(is.null(result$se_wcb))
  expect_s3_class(result$wcb_details, "lwdid_wcb_result")
  expect_identical(result$wcb_details$requested_n_bootstrap, 39L)
  expect_identical(result$wcb_details$weight_type, "mammen")
  expect_identical(result$wcb_details$restricted_model, "intercept_only")
  expect_false(isTRUE(result$wcb_details$use_fwildclusterboot))
})

test_that("TC-9.5.2: auto_wcb triggers on few clusters and emits a message", {
  panel <- make_e9_05_common_panel(n_units = 10L, n_periods = 8L, n_treated = 5L, n_clusters = 5L)
  result <- NULL
  captured_messages <- character(0)

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(lwdid(
      panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = "x1",
      vce = "cluster",
      cluster_var = "state",
      auto_wcb = TRUE,
      wcb_reps = 29L,
      use_fwildclusterboot = FALSE
    )),
    message = function(m) {
      captured_messages <<- c(captured_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_true(any(grepl("Wild Cluster Bootstrap", captured_messages, fixed = TRUE)))
  expect_identical(result$vce_type, "wild_cluster_bootstrap")
  expect_false(is.null(result$pvalue_wcb))
  expect_true(isTRUE(result$wcb_auto_triggered))
})

test_that("TC-9.5.6: quarter acts as a backward-compatible alias for season_var", {
  panel <- make_e9_05_seasonal_panel()

  result_alias <- local_suppress_lwdid_warnings(lwdid(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    d = "d",
    post = "post",
    rolling = "demeanq",
    quarter = "quarter",
    Q = 4L
  ))

  result_direct <- local_suppress_lwdid_warnings(lwdid(
    panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    d = "d",
    post = "post",
    rolling = "demeanq",
    season_var = "quarter",
    Q = 4L
  ))

  expect_s3_class(result_alias, "lwdid_result")
  expect_equal(result_alias$att, result_direct$att, tolerance = 1e-10)
  expect_equal(result_alias$se_att, result_direct$se_att, tolerance = 1e-10)
})

test_that("TC-9.5.3/5/7/10: seasonal integration preserves transformed inputs and clear errors", {
  seasonal_panel <- make_e9_05_seasonal_panel()

  seasonal_result <- local_suppress_lwdid_warnings(lwdid(
    seasonal_panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    d = "d",
    post = "post",
    rolling = "demeanq",
    season_var = "quarter",
    Q = 4L
  ))
  exclude_zero_result <- local_suppress_lwdid_warnings(lwdid(
    seasonal_panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    d = "d",
    post = "post",
    rolling = "demeanq",
    season_var = "quarter",
    Q = 4L,
    exclude_pre_periods = 0L
  ))
  exclude_two_result <- local_suppress_lwdid_warnings(lwdid(
    seasonal_panel,
    y = "y",
    ivar = "id",
    tvar = "time",
    d = "d",
    post = "post",
    rolling = "demeanq",
    season_var = "quarter",
    Q = 4L,
    exclude_pre_periods = 2L
  ))
  season_error <- tryCatch(
    local_suppress_lwdid_warnings(lwdid(
      seasonal_panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "d",
      post = "post",
      rolling = "demeanq"
    )),
    error = function(e) conditionMessage(e)
  )

  expect_identical(seasonal_result$rolling, "demeanq")
  expect_true(is.numeric(seasonal_result$att))
  expect_true(is.finite(unname(as.numeric(seasonal_result$att))))
  expect_identical(exclude_zero_result$exclude_pre_periods, 0L)
  expect_identical(exclude_two_result$exclude_pre_periods, 2L)
  expect_gt(
    abs(unname(as.numeric(exclude_zero_result$att)) - unname(as.numeric(exclude_two_result$att))),
    0
  )
  expect_match(season_error, "season_var", fixed = TRUE)
})

test_that("TC-9.5.3/9/11/14: package-facing WCB boundary, fallback message-only, and alpha CI stay stable", {
  common_panel <- make_e9_05_common_panel()

  explicit_result <- local_suppress_lwdid_warnings(lwdid(
    common_panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = "x1",
    vce = "bootstrap",
    cluster_var = "state",
    wcb_reps = 39L,
    wcb_seed = 17L,
    use_fwildclusterboot = FALSE
  ))
  alpha_95_result <- local_suppress_lwdid_warnings(lwdid(
    common_panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = "x1",
    vce = "bootstrap",
    cluster_var = "state",
    wcb_reps = 39L,
    wcb_seed = 17L,
    use_fwildclusterboot = FALSE,
    alpha = 0.05
  ))
  alpha_90_result <- local_suppress_lwdid_warnings(lwdid(
    common_panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = "x1",
    vce = "bootstrap",
    cluster_var = "state",
    wcb_reps = 39L,
    wcb_seed = 17L,
    use_fwildclusterboot = FALSE,
    alpha = 0.10
  ))
  fallback_capture <- capture_conditions(local_suppress_lwdid_warnings(lwdid(
    common_panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = "x1",
    vce = "bootstrap",
    cluster_var = "state",
    wcb_reps = 19L,
    wcb_seed = 7L,
    use_fwildclusterboot = TRUE
  )))
  cluster_error <- tryCatch(
    local_suppress_lwdid_warnings(lwdid(
      common_panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = "x1",
      vce = "bootstrap"
    )),
    error = function(e) conditionMessage(e)
  )

  expect_identical(explicit_result$.wcb_y_transformed, "y_trans_summary")
  expect_identical(nrow(explicit_result$.wcb_data), 12L)
  expect_true(all(explicit_result$.wcb_data$firstpost))
  expect_true(isTRUE(all.equal(alpha_95_result$att, alpha_90_result$att, tolerance = 1e-12)))
  expect_lt(
    unname(as.numeric(alpha_90_result$ci_upper_wcb - alpha_90_result$ci_lower_wcb)),
    unname(as.numeric(alpha_95_result$ci_upper_wcb - alpha_95_result$ci_lower_wcb))
  )
  expect_true(any(grepl("fwildclusterboot", fallback_capture$messages, fixed = TRUE)))
  expect_true(any(grepl("install.packages", fallback_capture$messages, fixed = TRUE)))
  expect_length(fallback_capture$warnings, 0L)
  expect_identical(fallback_capture$value$wcb_details$method, "native")
  expect_false(isTRUE(fallback_capture$value$wcb_details$use_fwildclusterboot))
  expect_match(cluster_error, "cluster_var", fixed = TRUE)
})

test_that("TC-9.5.4/12/13: WCB replaces top-level inference and print-summary expose WCB details", {
  panel <- make_e9_05_common_panel()

  result <- local_suppress_lwdid_warnings(lwdid(
    panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = "x1",
    vce = "bootstrap",
    cluster_var = "state",
    wcb_reps = 31L,
    wcb_type = "webb",
    wcb_seed = 77L,
    wcb_restricted_model = "with_controls",
    use_fwildclusterboot = FALSE
  ))

  expect_equal(result$se_att, result$se_wcb, tolerance = 1e-12)
  expect_equal(result$pvalue, result$pvalue_wcb, tolerance = 1e-12)
  expect_equal(result$ci_lower, result$ci_lower_wcb, tolerance = 1e-12)
  expect_equal(result$ci_upper, result$ci_upper_wcb, tolerance = 1e-12)

  print_output <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(print_output, "Wild Cluster Bootstrap")
  expect_match(print_output, "Weight type", fixed = TRUE)
  expect_match(print_output, "Requested bootstrap", fixed = TRUE)
  expect_match(print_output, "Actual bootstrap", fixed = TRUE)

  summary_output <- paste(capture.output(summary(result)), collapse = "\n")
  expect_match(summary_output, "Wild Cluster Bootstrap")
  expect_match(summary_output, "Weight type", fixed = TRUE)
  expect_match(summary_output, "Restricted model", fixed = TRUE)
})

test_that("TC-9.5.16a: explicit bootstrap with G<2 returns degenerate NaN WCB inference", {
  panel <- make_e9_05_common_panel(
    n_units = 8L,
    n_periods = 8L,
    n_treated = 4L,
    n_clusters = 1L
  )
  result <- NULL

  expect_no_error(
    result <- local_suppress_lwdid_warnings(lwdid(
      panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = "x1",
      vce = "bootstrap",
      cluster_var = "state",
      wcb_reps = 19L,
      use_fwildclusterboot = FALSE
    ))
  )

  expect_identical(result$vce_type, "wild_cluster_bootstrap")
  expect_true(is.nan(result$pvalue_wcb))
  expect_true(is.nan(result$ci_lower_wcb))
  expect_true(is.nan(result$ci_upper_wcb))
  expect_true(is.nan(result$se_wcb))
  expect_identical(result$wcb_details$n_clusters, 1L)
})

test_that("TC-9.5.16b: auto_wcb with G<2 also degrades to NaN WCB inference", {
  panel <- make_e9_05_common_panel(
    n_units = 8L,
    n_periods = 8L,
    n_treated = 4L,
    n_clusters = 1L
  )
  result <- NULL
  captured_messages <- character(0)

  expect_no_error(
    result <- withCallingHandlers(
      local_suppress_lwdid_warnings(lwdid(
        panel,
        y = "y",
        ivar = "id",
        tvar = "year",
        d = "d",
        post = "post",
        controls = "x1",
        vce = "cluster",
        cluster_var = "state",
        auto_wcb = TRUE,
        wcb_reps = 19L,
        use_fwildclusterboot = FALSE
      )),
      message = function(m) {
        captured_messages <<- c(captured_messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
  )

  expect_true(any(grepl("Wild Cluster Bootstrap", captured_messages, fixed = TRUE)))
  expect_true(isTRUE(result$wcb_auto_triggered))
  expect_true(is.nan(result$pvalue_wcb))
  expect_true(is.nan(result$se_wcb))
  expect_identical(result$wcb_details$n_clusters, 1L)
})

test_that("TC-9.5.15/17: non-triggering cluster paths expose stable FALSE auto-WCB flag", {
  few_cluster_panel <- make_e9_05_common_panel(
    n_units = 10L,
    n_periods = 8L,
    n_treated = 5L,
    n_clusters = 5L,
    seed = 99L
  )
  many_cluster_panel <- make_e9_05_common_panel(
    n_units = 40L,
    n_periods = 8L,
    n_treated = 20L,
    n_clusters = 20L,
    seed = 314L
  )
  few_cluster_messages <- character(0)
  many_cluster_messages <- character(0)

  few_cluster_result <- withCallingHandlers(
    local_suppress_lwdid_warnings(lwdid(
      few_cluster_panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = "x1",
      vce = "cluster",
      cluster_var = "state",
      auto_wcb = FALSE,
      use_fwildclusterboot = FALSE
    )),
    message = function(m) {
      few_cluster_messages <<- c(few_cluster_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  many_cluster_result <- withCallingHandlers(
    local_suppress_lwdid_warnings(lwdid(
      many_cluster_panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = "x1",
      vce = "cluster",
      cluster_var = "state",
      auto_wcb = TRUE,
      use_fwildclusterboot = FALSE
    )),
    message = function(m) {
      many_cluster_messages <<- c(many_cluster_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_identical(few_cluster_result$vce_type, "cluster")
  expect_identical(few_cluster_result$wcb_auto_triggered, FALSE)
  expect_null(few_cluster_result$wcb_details)
  expect_length(few_cluster_messages, 0L)

  expect_identical(many_cluster_result$vce_type, "cluster")
  expect_identical(many_cluster_result$wcb_auto_triggered, FALSE)
  expect_null(many_cluster_result$wcb_details)
  expect_length(many_cluster_messages, 0L)
})
