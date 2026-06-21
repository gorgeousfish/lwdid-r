# lwdid 0.1.0

## Initial release candidate

* Added Lee-Wooldridge difference-in-differences estimators for common timing
  and staggered adoption designs.
* Added regression adjustment, inverse probability weighting, doubly robust
  IPWRA, and propensity score matching estimators.
* Added robust variance options, including heteroskedasticity-consistent and
  cluster-robust covariance estimators.
* Added randomization inference and a native wild cluster bootstrap path, with
  optional integration for `fwildclusterboot` when available.
* Added diagnostics for pre-treatment trends, overlap, clustering, and
  selection, together with print, summary, plotting, and export helpers.
* Added package vignettes, examples, and test coverage used by the accompanying
  R Journal manuscript replication materials.
