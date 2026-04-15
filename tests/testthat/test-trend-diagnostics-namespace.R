library(testthat)

test_that("E8-06.4: NAMESPACE exports trend diagnostics public surface", {
  namespace_lines <- readLines(resolve_package_source_file("NAMESPACE"), warn = FALSE)

  expected_entries <- c(
    "export(lwdid_test_parallel_trends)",
    "export(lwdid_diagnose_heterogeneous_trends)",
    "export(lwdid_recommend_transformation)",
    "export(plot_cohort_trends)",
    "S3method(print,lwdid_parallel_trends)",
    "S3method(summary,lwdid_parallel_trends)",
    "S3method(plot,lwdid_parallel_trends)",
    "S3method(print,lwdid_heterogeneous_trends)",
    "S3method(summary,lwdid_heterogeneous_trends)",
    "S3method(plot,lwdid_heterogeneous_trends)",
    "S3method(print,lwdid_transformation_recommendation)",
    "S3method(summary,lwdid_transformation_recommendation)",
    "S3method(plot,lwdid_transformation_recommendation)"
  )

  missing_entries <- setdiff(expected_entries, namespace_lines)

  expect(
    length(missing_entries) == 0L,
    paste(
      "Missing NAMESPACE entries:",
      paste(missing_entries, collapse = ", ")
    )
  )
})
