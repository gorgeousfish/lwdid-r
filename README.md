# lwdid

**Lee-Wooldridge Difference-in-Differences Estimation for R**

[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE.md)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)](NEWS.md)

## Overview

`lwdid` implements the **Rolling Difference-in-Differences Estimator** proposed by Lee and Wooldridge ([2025](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4516518), [2026](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686)) for R. Through unit-specific time-series transformations (demeaning or detrending), the package converts panel DiD estimation into standard cross-sectional treatment effect problems, enabling several treatment-effect estimators and the package's supported inference options.

The package supports common timing and staggered adoption treatment settings in panel data. It keeps estimation, diagnostics, plots, exports, and replication-oriented metadata attached to fitted R objects so users can inspect the analysis path rather than copy results across disconnected scripts.

**Features**:

- **Two designs**: Common Timing (simultaneous treatment) and Staggered Adoption (cohort-specific timing)
- **Four estimators**: RA (Regression Adjustment), IPW (Inverse Probability Weighting), IPWRA (Doubly Robust), PSM (Propensity Score Matching)
- **Inference options**: Homoskedastic / HC0-HC4 / cluster-robust standard errors, Fisher Randomization Inference, Wild Cluster Bootstrap
- **Rich diagnostics**: Parallel trends tests, clustering diagnostics, selection mechanism diagnostics, sensitivity analysis
- **Data preprocessing**: Individual-level data aggregation to panel format via `aggregate_to_panel()`
- **Implementation**: `data.table`-oriented computation with optional parallel processing via `future`
- **Export**: LaTeX tables (`to_latex()`, `to_latex_comparison()`), CSV (`to_csv()`), event study plots

## Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Transformations](#transformations)
- [Estimators](#estimators)
- [Aggregation](#aggregation)
- [Variance-Covariance Estimation](#variance-covariance-estimation)
- [Diagnostics](#diagnostics)
- [Inference](#inference)
- [Post-Estimation Tools](#post-estimation-tools)
- [Built-in Datasets](#built-in-datasets)
- [Citation](#citation)
- [Authors](#authors)
- [See Also](#see-also)

## Installation

**Requirements:** R (>= 4.0.0)

### From GitHub (development version)

```r
# install.packages("devtools")
devtools::install_github("gorgeousfish/lwdid-r")
```

### From CRAN (planned)

```r
install.packages("lwdid")
```

## Optional Dependencies

The core functionality requires only `data.table` and `sandwich`.
Optional packages unlock additional features:

```r
# Wild Cluster Bootstrap inference (archived from CRAN; install from R-universe)
install.packages("fwildclusterboot",
                 repos = "https://s3alfisc.r-universe.dev")

# Propensity Score Matching (PSM estimator)
install.packages(c("MatchIt", "cobalt"))

# Parallel computation (RI / Bootstrap / Staggered loops)
install.packages(c("future", "future.apply"))

# Visualization
install.packages("ggplot2")
```

## Quick Start

### Example 1: Common Timing (Smoking Data)

The `smoking` dataset contains the California Proposition 99 tobacco control study (1970–2000). California is the sole treated state.

```r
library(lwdid)
data(smoking)

result <- lwdid(
  data = smoking,
  y = "lcigsale",
  ivar = "state",
  tvar = "year",
  d = "d",
  post = "post",
  rolling = "demean"
)

summary(result)
```

### Example 2: Staggered Adoption (Castle Doctrine Data)

The `castle` dataset examines Castle Doctrine laws and homicide rates. States adopted the law at different times between 2005 and 2009.

```r
data(castle)

result <- lwdid(
  data = castle,
  y = "lhomicide",
  ivar = "sid",
  tvar = "year",
  gvar = "gvar",
  rolling = "demean",
  aggregate = "overall"
)

summary(result)
```

For `aggregate = "overall"`, the reported ATT is a cohort-size weighted
average of the cohort-time effects. The result object and exports retain
`cohort_weights`; `effective_weights` record the weights after any rows dropped
from the regression sample.

### Example 3: Event Study

```r
result_es <- lwdid(
  data = castle,
  y = "lhomicide",
  ivar = "sid",
  tvar = "year",
  gvar = "gvar",
  rolling = "demean",
  aggregate = "event_time"
)

summary(result_es)
plot(result_es)
```

For `aggregate = "event_time"`, the printed scalar ATT summarizes valid
cohort-period `(g, r)` effects. The event-study table and plot report WATT(e)
rows by relative time, with cohort-size weights renormalized over the cohorts
available at each event time. WATT(e) standard errors use the stored
cohort-specific standard errors with diagonal aggregation; cross-cohort
covariance is not modeled unless a future joint covariance representation is
added.

Use `extract_effects(result_es, type = "event_time_contributions")` to inspect
the cohort weights and cohort-specific contributions behind each WATT(e) row.
The same contribution table can be written with
`to_csv(result_es, file = "event-time-contributions.csv", what = "event_time_contributions")`
or read from `to_dict(result_es)$event_time_contributions`.
Use `to_latex(result_es, include_event_time = TRUE)` to export the WATT(e)
rows with standard errors, confidence intervals, contributing-cohort counts,
and compact support/overlap metadata.
Plot data returned by `plot_event_study(result_es, return_data = TRUE)`
carries the same standard-error metadata for the plotted event-time rows,
including `se_aggregation` and `covariance_assumption`. For non-RA
event-time results, the returned plot data also preserves available overlap
and support-count summaries, such as propensity-score ranges, weight
dispersion, and treated/control counts.

```r
# Hide pre-treatment periods; the default -1 anchor is omitted as well.
plot(result_es, show_pre_treatment = FALSE)

# Normalize to an observed event time; this does not relabel the observed reference row as a visual anchor.
plot(result_es, ref_period = 0)
```

When `ref_period` is used, returned plot data suppresses pointwise p-values
and significance flags. The plotted values are shifted contrasts relative to
the reference period, and those contrasts require event-time contrast
covariance rather than the original pointwise standard errors.

### Example 4: Using Controls and IPW Estimator

```r
result_ipw <- lwdid(
  data = castle,
  y = "lhomicide",
  ivar = "sid",
  tvar = "year",
  gvar = "gvar",
  rolling = "demean",
  estimator = "ipw",
  ps_controls = c("police", "unemployrt", "income"),
  aggregate = "overall"
)

summary(result_ipw)
```

## Transformations

| `rolling =`  | Description                          | Min. pre-treatment periods |
| -------------- | ------------------------------------ | -------------------------- |
| `"demean"`   | Subtract pre-treatment mean          | 1                          |
| `"detrend"`  | Remove unit-specific linear trend    | 2                          |
| `"demeanq"`  | Seasonal demeaning (quarterly data)  | 1 per season               |
| `"detrendq"` | Seasonal detrending (quarterly data) | 2 per season               |

## Estimators

| `estimator =` | Description                   | Controls required | Doubly robust |
| --------------- | ----------------------------- | :---------------: | :-----------: |
| `"ra"`        | Regression Adjustment         |     Optional     |      No      |
| `"ipw"`       | Inverse Probability Weighting |        Yes        |      No      |
| `"ipwra"`     | Doubly Robust (IPW + RA)      |        Yes        |      Yes      |
| `"psm"`       | Propensity Score Matching     |        Yes        |      No      |

Use `controls` for RA outcome-regression controls. For IPW and PSM, prefer
`ps_controls` to state the propensity-score specification directly; if
`ps_controls` is `NULL`, they fall back to `controls`. IPWRA uses `controls`
for the outcome model and `ps_controls` for the propensity-score model, again
the propensity model will fall back to `controls` when `ps_controls` is
omitted.

## Aggregation

For staggered adoption designs, treatment effects are first estimated for each cohort-time pair and then aggregated:

| `aggregate =`     | Description                                                |
| ------------------- | ---------------------------------------------------------- |
| `"cohort"`        | Group-time ATT by treatment cohort (default)               |
| `"overall"`       | Single ATT weighted by treated-unit cohort sizes           |
| `"event_time"`    | WATT(e) by relative time; scalar ATT remains a `(g,r)` summary |
| `"calendar_time"` | ATT by calendar time period                                |
| `"att_gt"`        | Return all individual group-time ATT(g,t) estimates        |

The `control_group` parameter controls which untreated units serve as comparisons:

| `control_group =`   | Description                                                           |
| --------------------- | --------------------------------------------------------------------- |
| `"not_yet_treated"` | Units not yet treated at time*t* (default; larger comparison group) |
| `"never_treated"`   | Only units that are never treated during the sample period            |

## Variance-Covariance Estimation

| `vce =`                | Description                                |
| ------------------------ | ------------------------------------------ |
| `NULL`                 | OLS (homoskedastic)                        |
| `"robust"` / `"hc1"` | HC1 heteroskedasticity-robust              |
| `"hc0"` – `"hc4"`   | HC-family (HC3 recommended for small*N*) |
| `"cluster"`            | Cluster-robust (requires `cluster_var`)  |

## Diagnostics

```r
# Parallel trends test (staggered)
pt <- lwdid_test_parallel_trends(
  data = castle, y = "lhomicide",
  ivar = "sid", tvar = "year", gvar = "gvar"
)
summary(pt)

# Sensitivity analysis (common timing)
sens <- lwdid_sensitivity(
  data = smoking, y = "lcigsale",
  ivar = "state", tvar = "year",
  d = "d", post = "post",
  type = "all"
)
summary(sens)

# Full diagnostic suite (staggered)
diag <- lwdid_diagnose(
  data = castle, y = "lhomicide",
  ivar = "sid", tvar = "year", gvar = "gvar"
)
summary(diag)
```

## Inference

### Randomization Inference

Randomization inference (RI) provides randomization-based p-values without relying on asymptotic normality, which is particularly useful in small-*N* settings:

```r
result_ri <- lwdid(
  data = smoking,
  y = "lcigsale",
  ivar = "state",
  tvar = "year",
  d = "d",
  post = "post",
  rolling = "demean",
  ri = TRUE,
  ri_method = "permutation",
  rireps = 1000,
  seed = 42
)

summary(result_ri)
```

RI output reports both `ri_pvalue` and `ri_observed_stat`. The p-value is
computed against `ri_observed_stat`, not necessarily the scalar `att` printed at
the top of a staggered result. For common-timing designs these coincide. For
staggered designs, inspect `ri_target`: `aggregate = "overall"` targets the
overall ATT, `aggregate = "cohort"` targets the selected cohort ATT, and
`aggregate = "none"` targets the first finite cohort-time ATT.

### Wild Cluster Bootstrap

Use `vce = "bootstrap"` to request Wild Cluster Bootstrap (WCB) inference directly. When `vce = "cluster"` is used, `lwdid` automatically upgrades to WCB only for fewer than 20 clusters.

```r
result_wcb <- lwdid(
  data = smoking,
  y = "lcigsale",
  ivar = "state",
  tvar = "year",
  d = "d",
  post = "post",
  rolling = "demean",
  vce = "bootstrap",
  cluster_var = "state",
  wcb_reps = 999,
  wcb_seed = 42,
  use_fwildclusterboot = FALSE
)

summary(result_wcb)
```

WCB output reports both requested and actual bootstrap draws. With
Rademacher weights and at most 12 clusters, `lwdid` uses full sign-pattern
enumeration, so the actual draw count can differ from `wcb_reps`.

## Post-Estimation Tools

### Multi-Specification Comparison

Compare results across different estimators, transformations, or control group strategies:

```r
# Estimate with different specifications
res_ra    <- lwdid(data, y, ivar, tvar, gvar, estimator = "ra")
res_ipwra <- lwdid(data, y, ivar, tvar, gvar, estimator = "ipwra")
res_ipw   <- lwdid(data, y, ivar, tvar, gvar, estimator = "ipw")

# Side-by-side comparison
compare(RA = res_ra, IPWRA = res_ipwra, IPW = res_ipw)
```

### Integration with modelsummary and broom

lwdid results work seamlessly with the tidyverse ecosystem:

```r
library(modelsummary)

# Publication-ready tables
modelsummary(list(RA = res_ra, IPWRA = res_ipwra))

# Tidy output for custom processing
tidy(result)                    # Coefficient estimates
tidy(result, type = "effects")  # Period-by-period effects
glance(result)                  # Model-level statistics
```

### fixest-Style Accessors

For users familiar with fixest workflows:

```r
se(result)        # Standard errors
tstat(result)     # t-statistics
pvalue(result)    # p-values
coeftable(result) # Full coefficient table
```

## Built-in Datasets

| Dataset     | Design        | Description                                                 | Obs    | Units          |
| ----------- | ------------- | ----------------------------------------------------------- | ------ | -------------- |
| `smoking` | Common Timing | California Proposition 99 tobacco control (1970–2000)      | 1,209  | 39 states      |
| `castle`  | Staggered     | Castle Doctrine / Stand Your Ground laws (2000–2010)       | 550    | 50 states      |
| `walmart` | Staggered     | Walmart store openings and county labor markets (1977–1999) | 29,440 | 1,280 counties |

## Citation

If you use `lwdid` in your research, please cite the following papers:

**APA Format:**

> Cai, X., & Xu, W. (2026). *lwdid: Lee-Wooldridge Difference-in-Differences Estimation for R* (Version 0.1.0) [Computer software]. GitHub. https://github.com/gorgeousfish/lwdid-r
>
> Lee, S. J., & Wooldridge, J. M. (2025). A simple transformation approach to difference-in-differences estimation for panel data. Available at SSRN 4516518.
>
> Lee, S. J., & Wooldridge, J. M. (2026). Simple Approaches to Inference with Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes. Available at SSRN 5325686.

**BibTeX:**

```bibtex
@software{lwdid2026r,
  title={lwdid: Lee-Wooldridge Difference-in-Differences Estimation for R},
  author={Xuanyu Cai and Wenli Xu},
  year={2026},
  version={0.1.0},
  url={https://github.com/gorgeousfish/lwdid-r}
}

@article{lee2025simple,
  title={A Simple Transformation Approach to Difference-in-Differences
         Estimation for Panel Data},
  author={Lee, Soo Jeong and Wooldridge, Jeffrey M.},
  year={2025},
  note={Available at SSRN 4516518}
}

@article{lee2026simple,
  title={Simple Approaches to Inference with Difference-in-Differences
         Estimators with Small Cross-Sectional Sample Sizes},
  author={Lee, Soo Jeong and Wooldridge, Jeffrey M.},
  year={2026},
  note={Available at SSRN 5325686}
}
```

## Authors

**R Implementation:**

- **Xuanyu Cai**, City University of Macau
  Email: [xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
- **Wenli Xu**, City University of Macau
  Email: [wlxu@cityu.edu.mo](mailto:wlxu@cityu.edu.mo)

**Methodology:**

- **[Soo Jeong Lee](https://sites.google.com/view/sjlee-econ/home)**, Michigan State University
- **Jeffrey M. Wooldridge**, Michigan State University

## License

AGPL-3.0. See [LICENSE.md](LICENSE.md) for details.

## See Also

- Original Stata package by Lee and Wooldridge: [https://github.com/Soo-econ/lwdid](https://github.com/Soo-econ/lwdid)
- Python implementation: [https://github.com/gorgeousfish/lwdid-py](https://github.com/gorgeousfish/lwdid-py)
