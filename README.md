# lwdid

**Lee-Wooldridge Difference-in-Differences Estimation for R**

[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE.md)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)]()

![lwdid](image/image.png)

## Overview

`lwdid` implements the **Rolling Difference-in-Differences Estimator** proposed by Lee and Wooldridge ([2025](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4516518), [2026](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686)) for R. Through unit-specific time-series transformations (demeaning or detrending), the package converts panel DiD estimation into standard cross-sectional treatment effect problems, enabling application of multiple estimators and exact small-sample inference.

The package provides fast and flexible estimation for both common timing and staggered adoption treatment settings in panel data, and covers standard large-*N* asymptotic inference as well as small-*N* settings where conventional inference may not be reliable.

**Features**:

- **Two designs**: Common Timing (simultaneous treatment) and Staggered Adoption (cohort-specific timing)
- **Four estimators**: RA (Regression Adjustment), IPW (Inverse Probability Weighting), IPWRA (Doubly Robust), PSM (Propensity Score Matching)
- **Complete inference**: Homoskedastic / HC0–HC4 / cluster-robust standard errors, Fisher Randomization Inference, Wild Cluster Bootstrap
- **Rich diagnostics**: Parallel trends tests, clustering diagnostics, selection mechanism diagnostics, sensitivity analysis
- **Data preprocessing**: Individual-level data aggregation to panel format via `aggregate_to_panel()`
- **High performance**: Vectorized `data.table`-based computation with optional parallel processing via `future`
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
- [Built-in Datasets](#built-in-datasets)
- [Vignettes](#vignettes)
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

plot(result_es)
```

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
  controls = c("police", "unemployrt", "income"),
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

## Aggregation

For staggered adoption designs, treatment effects are first estimated for each cohort-time pair and then aggregated:

| `aggregate =`     | Description                                                |
| ------------------- | ---------------------------------------------------------- |
| `"cohort"`        | Group-time ATT by treatment cohort (default)               |
| `"overall"`       | Single weighted average ATT across all cohorts and periods |
| `"event_time"`    | ATT by relative time to treatment (for event study plots)  |
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

Randomization inference (RI) provides exact p-values without relying on asymptotic normality, which is particularly useful in small-*N* settings:

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
  rireps = 1000,
  seed = 42
)

summary(result_ri)
```

### Wild Cluster Bootstrap

Wild Cluster Bootstrap (WCB) is automatically applied when `vce = "cluster"` is specified, providing asymptotic refinement for cluster-robust inference:

```r
result_wcb <- lwdid(
  data = castle,
  y = "lhomicide",
  ivar = "sid",
  tvar = "year",
  gvar = "gvar",
  rolling = "demean",
  vce = "cluster",
  cluster_var = "sid",
  aggregate = "overall"
)

summary(result_wcb)
```

## Built-in Datasets

| Dataset     | Design        | Description                                                 | Obs    | Units          |
| ----------- | ------------- | ----------------------------------------------------------- | ------ | -------------- |
| `smoking` | Common Timing | California Proposition 99 tobacco control (1970–2000)      | 1,209  | 39 states      |
| `castle`  | Staggered     | Castle Doctrine / Stand Your Ground laws (2000–2010)       | 550    | 50 states      |
| `walmart` | Staggered     | Walmart store openings and local labor markets (1977–1999) | 29,440 | 1,280 counties |

## Vignettes

- `vignette("getting-started")`: Data format, transformations, common timing and staggered examples
- `vignette("advanced-usage")`: Multiple estimators, sensitivity, diagnostics, WCB, RI, parallel computation, export
- `vignette("methodology")`: Mathematical framework, identification, and comparison with Callaway & Sant'Anna (2021)

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
