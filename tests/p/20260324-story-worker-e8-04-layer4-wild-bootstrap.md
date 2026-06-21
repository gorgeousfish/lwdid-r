# E8-04 Layer 4 wild-bootstrap convergence

## 本轮结论

- `wild_bootstrap_coverage_small_g` 当前已具备 executable regression：
  `coverage_rate = 0.96`，满足 contract 的 `> 0.85` acceptance band。
- `wild_bootstrap_relative_gap` 当前已具备 executable regression：
  `coverage_std_rate = 0.91`、`coverage_wild_rate = 0.96`、
  `coverage_rate_gap = 0.05`，满足 contract 的 `>= -0.10` acceptance band。
- 两条 scenario 都已锁定 Python effective WCB contract，而不是“199 次随机
  bootstrap”：
  `requested_n_bootstrap = 199`、`actual_n_bootstrap = 1024`、
  `weight_type = rademacher`、`impose_null = TRUE`、
  `ci_method = percentile_t`、`full_enumeration = TRUE`。

## 关键证据

- 新增 RED/GREEN 回归：
  - `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`
- 新增最小 WCB helper：
  - `lwdid-r/R/clustering_diagnostics.R`
- fresh rerun：
  - `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics-parity", reporter="summary")'`
  - `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics", reporter="summary")'`

## 影响

- `story-E8-04` 的 Layer 4 现已完成全部 9 个 Monte Carlo scenario 的
  executable parity。
- 当前 parity blocker 已清零；后续若仍保留 `story-E8-04` 为 active story，
  应把焦点转到 fresh `devtools::check()` 与 controller closure audit，
  不再回退到 wild-bootstrap pending。
