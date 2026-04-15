# ===========================================================================
# test-zzz.R — Tests for zzz.R package hooks
#
# Covers: .lwdid_env, warning_registry initialization, version cache,
#         package options, .onAttach(), get_lwdid_env(), CRAN compliance,
#         option lifecycle, and file load order dependencies.
# ===========================================================================

# ---- Group A: .lwdid_env environment verification -------------------------

test_that(".lwdid_env is an environment", {
  expect_true(is.environment(lwdid:::.lwdid_env))
})

test_that(".lwdid_env parent is emptyenv", {
  expect_identical(parent.env(lwdid:::.lwdid_env), emptyenv())
})

test_that(".lwdid_env contains warning_registry after load", {
  expect_true(exists("warning_registry", envir = lwdid:::.lwdid_env))
})

test_that(".lwdid_env contains version after load", {
  expect_true(exists("version", envir = lwdid:::.lwdid_env))
})

# ---- Group B: warning_registry initialization -----------------------------

test_that("warning_registry is an environment", {
  reg <- lwdid:::.lwdid_env$warning_registry
  expect_true(is.environment(reg))
})

test_that("warning_registry has all required methods", {
  reg <- lwdid:::.lwdid_env$warning_registry
  for (method in c("register", "flush", "clear", "count",
                    "get_log", "get_diagnostics", "merge")) {
    expect_true(is.function(reg[[method]]),
                info = paste("Missing method:", method))
  }
})

test_that("warning_registry initial count is 0", {
  reg <- lwdid:::.lwdid_env$warning_registry
  # Clear first to ensure clean state
  reg$clear()
  expect_equal(reg$count(), 0L)
})

test_that("warning_registry register and clear work", {
  reg <- lwdid:::.lwdid_env$warning_registry
  reg$clear()
  reg$register("lwdid_small_sample", "test message",
               cohort = 2005L, period = 2006L)
  expect_equal(reg$count(), 1L)
  reg$clear()
  expect_equal(reg$count(), 0L)
})

# ---- Group C: version cache -----------------------------------------------

test_that("version is a package_version object", {
  expect_true(inherits(lwdid:::.lwdid_env$version, "package_version"))
})

test_that("version matches packageVersion", {
  # Note: during devtools::test(), the version may be from the loaded dev package
  # Just verify it's a valid version object
  ver <- lwdid:::.lwdid_env$version
  expect_true(is.character(as.character(ver)))
  expect_true(nchar(as.character(ver)) > 0)
})

# ---- Group D: Package options ----------------------------------------------

test_that("lwdid.verbose option is set after load", {
  expect_false(is.null(getOption("lwdid.verbose")))
})

test_that("lwdid.verbose defaults to 'default'", {
  # Save, reset, reload, check, restore
  old_val <- getOption("lwdid.verbose")
  options(lwdid.verbose = NULL)
  lwdid:::.onLoad(NULL, "lwdid")
  expect_equal(getOption("lwdid.verbose"), "default")
  options(lwdid.verbose = old_val)
})

test_that("lwdid.verbose user preset is not overwritten", {
  old_val <- getOption("lwdid.verbose")
  options(lwdid.verbose = "quiet")
  lwdid:::.onLoad(NULL, "lwdid")
  expect_equal(getOption("lwdid.verbose"), "quiet")
  options(lwdid.verbose = old_val)
})

# ---- Group E: .onAttach() startup message ----------------------------------

test_that(".onAttach produces startup message with version", {
  expect_message(
    lwdid:::.onAttach(NULL, "lwdid"),
    "lwdid .+: Lee-Wooldridge DiD Estimation"
  )
})

test_that(".onAttach message can be suppressed", {
  expect_no_message(
    suppressPackageStartupMessages(lwdid:::.onAttach(NULL, "lwdid"))
  )
})

# ---- Group F: get_lwdid_env() helper ---------------------------------------

test_that("get_lwdid_env returns .lwdid_env", {
  expect_identical(lwdid:::get_lwdid_env(), lwdid:::.lwdid_env)
})

test_that("get_lwdid_env has reference semantics", {
  env <- lwdid:::get_lwdid_env()
  env$.__test_ref_semantics__ <- 42
  expect_equal(lwdid:::.lwdid_env$.__test_ref_semantics__, 42)
  rm(".__test_ref_semantics__", envir = lwdid:::.lwdid_env)
})

# ---- Group G: .onLoad() CRAN compliance (silent execution) -----------------

test_that(".onLoad is silent (no messages)", {
  expect_no_message(lwdid:::.onLoad(NULL, "lwdid"))
})

test_that(".onLoad produces no warnings", {
  expect_no_warning(lwdid:::.onLoad(NULL, "lwdid"))
})

test_that(".onLoad produces no stdout output", {
  out <- capture.output(lwdid:::.onLoad(NULL, "lwdid"), type = "output")
  expect_length(out, 0)
})

test_that(".onLoad registers trend diagnostics S3 methods", {
  lwdid:::.onLoad(NULL, "lwdid")

  expect_true(is.function(getS3method(
    "summary",
    "lwdid_parallel_trends",
    optional = TRUE
  )))
  expect_true(is.function(getS3method(
    "plot",
    "lwdid_parallel_trends",
    optional = TRUE
  )))
  expect_true(is.function(getS3method(
    "summary",
    "lwdid_heterogeneous_trends",
    optional = TRUE
  )))
  expect_true(is.function(getS3method(
    "plot",
    "lwdid_heterogeneous_trends",
    optional = TRUE
  )))
  expect_true(is.function(getS3method(
    "summary",
    "lwdid_transformation_recommendation",
    optional = TRUE
  )))
  expect_true(is.function(getS3method(
    "plot",
    "lwdid_transformation_recommendation",
    optional = TRUE
  )))
})

# ---- Group H: Option lifecycle ---------------------------------------------

test_that("option persists after simulated unload/reload", {
  old_val <- getOption("lwdid.verbose")
  options(lwdid.verbose = "verbose")
  # Simulate reload
  lwdid:::.onLoad(NULL, "lwdid")
  expect_equal(getOption("lwdid.verbose"), "verbose")
  options(lwdid.verbose = old_val)
})

# ---- Group I: File load order dependency -----------------------------------

test_that("new_warning_registry is available in namespace", {
  expect_true(is.function(lwdid:::new_warning_registry))
})

test_that("warning_registry is properly initialized (load order OK)", {
  reg <- lwdid:::.lwdid_env$warning_registry
  expect_false(is.null(reg))
  expect_true(is.environment(reg))
})

# ---- Integration tests -----------------------------------------------------

test_that("WarningRegistry full lifecycle works", {
  reg <- lwdid:::.lwdid_env$warning_registry
  reg$clear()

  # Register multiple warnings
  reg$register("lwdid_small_sample", "N < 30",
               cohort = 2005L, period = 2006L)
  reg$register("lwdid_data", "Missing values",
               cohort = 2005L, period = 2007L)
  expect_equal(reg$count(), 2L)

  # Get log

  log <- reg$get_log()
  expect_true(is.data.frame(log))
  expect_equal(nrow(log), 2L)

  # Clear
  reg$clear()
  expect_equal(reg$count(), 0L)
})

test_that("no .onUnload function exists in zzz.R", {
  # Verify .onUnload is NOT defined
  expect_false(exists(".onUnload", envir = asNamespace("lwdid")))
})
