## R CMD check results

0 errors | 0 warnings | 2 notes

* NOTE: New submission
  - This is a first submission of the package.
  
* NOTE: License with restrictions
  - AGPL-3 + file LICENSE is used intentionally to ensure all derivative 
    works remain open source, which is important for reproducibility in 
    academic research. The LICENSE file contains the standard AGPL-3 text.

## Test environments

* local macOS Sequoia 15.3 (aarch64-apple-darwin), R 4.4.x
* GitHub Actions: ubuntu-latest (release), windows-latest (release), macOS-latest (release)

## Downstream dependencies

This is a new package submission. There are currently no downstream dependencies.

## Package purpose

This package implements the Lee and Wooldridge (2025, 2026) simple 
transformation approach for difference-in-differences (DiD) estimation 
with panel data. It provides a unified framework for both common timing 
and staggered adoption designs with multiple estimators.

References:
- Lee, S. and Wooldridge, J.M. (2025). "A Simple Transformation Approach 
  to Difference-in-Differences Estimation for Panel Data." Available at 
  SSRN: https://ssrn.com/abstract=4516518
- Lee, S. and Wooldridge, J.M. (2026). "Staggered DiD Designs with 
  Heterogeneous Trends." Available at SSRN: https://ssrn.com/abstract=5325686

## External URLs

All URLs in the package documentation have been verified:
- GitHub repository: https://github.com/gorgeousfish/lwdid-r
- Issue tracker: https://github.com/gorgeousfish/lwdid-r/issues
## R CMD check results

Local staged source-tarball check after release metadata and installed-size
cleanup:

* This is an initial release candidate.
* Local staged `R CMD check` with tests and examples enabled, and with manual
  and vignette rebuild checks disabled: 0 errors, 0 warnings, 0 notes.
* Lazy data are compressed with `LazyDataCompression: xz`; the staged
  installed-size check records `* checking installed package size ... OK`.
* The package no longer declares the optional non-CRAN `fwildclusterboot`
  backend in `Enhances` or `Additional_repositories`. The code path detects
  that backend dynamically when it is installed and otherwise uses the native
  wild-cluster-bootstrap implementation.
* External URLs should be rechecked close to release submission because the
  current local URL evidence contains network timeouts rather than confirmed
  broken links.
* A final CRAN submission still requires the final release environment and any
  external pre-submission checks chosen for the release.

## Test environments

* local macOS staged installed-size check: `Status: OK` under the same
  manual/vignette boundary as the staged release check.
* local macOS full staged source-tarball check: `Status: OK` with tests and
  examples enabled and manual/vignette rebuild checks disabled.

## Submission status

This file records the current release-candidate intent after dependency metadata
cleanup and staged source-tarball checking. A final CRAN submission still
requires release-time URL rechecks and any external CRAN pre-submission checks
chosen for the release.

## Downstream dependencies

There are currently no known downstream dependencies for this initial release
candidate.
