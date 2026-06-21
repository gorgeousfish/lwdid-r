# This file is part of the standard testthat setup.
# It is run during `R CMD check` to execute all tests.

library(testthat)
library(lwdid)

evidence_filter <- paste(
  c(
    "parity",
    "realdata",
    "monte-carlo",
    "layer4-contract",
    "comprehensive",
    "doc-contract",
    "seasonal-wcb",
    "wild-bootstrap",
    "lwdid-integration",
    "sensitivity-diagnostics",
    "aggregate",
    "demeanq",
    "non-ra-bootstrap-diagnostics",
    "plot-event-study",
    "plot-diagnostics",
    "results",
    "stata-consistency"
  ),
  collapse = "|"
)

if (identical(Sys.getenv("LWDID_RUN_EVIDENCE_TESTS"), "true")) {
  test_check("lwdid")
} else {
  test_check("lwdid", filter = evidence_filter, invert = TRUE)
}
