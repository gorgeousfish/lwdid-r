# ============================================================================
# test-wcb-export-contract.R
# WCB result exports must preserve the inference settings needed for reporting.
# ============================================================================

test_that("WCB exports preserve real-data bootstrap metadata", {
  data("smoking", package = "lwdid")

  result <- suppressWarnings(lwdid(
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
  ))

  summary_text <- paste(capture.output(summary(result)), collapse = "\n")
  expect_match(summary_text, "Wild Cluster Bootstrap", fixed = TRUE)
  expect_match(summary_text, "Requested bootstrap draws: 99", fixed = TRUE)
  expect_match(summary_text, "Actual bootstrap draws: 99", fixed = TRUE)

  latex <- to_latex(result, stars = FALSE)
  expect_match(latex, "WCB weight & rademacher", fixed = TRUE)
  expect_match(latex, "WCB clusters & 39", fixed = TRUE)
  expect_match(latex, "WCB draws & 99 actual / 99 requested", fixed = TRUE)
  expect_match(latex, "WCB restricted model & with\\_controls", fixed = TRUE)

  csv_file <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_file), add = TRUE)
  csv <- suppressMessages(to_csv(result, file = csv_file, what = "summary"))

  expect_equal(csv$pvalue_wcb, result$pvalue_wcb)
  expect_equal(csv$se_wcb, result$se_wcb)
  expect_equal(csv$wcb_weight_type, "rademacher")
  expect_equal(csv$wcb_n_clusters, 39L)
  expect_equal(csv$wcb_requested_n_bootstrap, 99L)
  expect_equal(csv$wcb_actual_n_bootstrap, 99L)
  expect_equal(csv$wcb_restricted_model, "with_controls")
  expect_false(csv$wcb_auto_triggered)

  dict <- to_dict(result)
  expect_equal(dict$pvalue_wcb, result$pvalue_wcb)
  expect_equal(dict$se_wcb, result$se_wcb)
  expect_equal(dict$wcb_weight_type, "rademacher")
  expect_equal(dict$wcb_n_clusters, 39L)
  expect_equal(dict$wcb_requested_n_bootstrap, 99L)
  expect_equal(dict$wcb_actual_n_bootstrap, 99L)
  expect_equal(dict$wcb_restricted_model, "with_controls")
  expect_false(dict$wcb_auto_triggered)
})

test_that("WCB comparison tables preserve per-model bootstrap metadata", {
  data("smoking", package = "lwdid")

  base_result <- suppressWarnings(lwdid(
    data = smoking,
    y = "lcigsale",
    ivar = "state",
    tvar = "year",
    d = "d",
    post = "post",
    rolling = "demean"
  ))

  wcb_result <- suppressWarnings(lwdid(
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
  ))

  comparison <- to_latex_comparison(
    base_result,
    wcb_result,
    model_names = c("Analytic", "WCB"),
    stars = FALSE
  )

  expect_match(comparison, "VCE & homoskedastic & Wild Cluster Bootstrap", fixed = TRUE)
  expect_match(comparison, "WCB weight & -- & rademacher", fixed = TRUE)
  expect_match(comparison, "WCB clusters & -- & 39", fixed = TRUE)
  expect_match(comparison, "WCB draws & -- & 99 actual / 99 requested", fixed = TRUE)
})
