# E8-04 Layer 4 coverage_large_g regression

## 本轮结论

- `confirmed`: `story-E8-04` 已不再是“Layer 4 完全没有 executable regression”。
  `coverage_large_g` 已在 R 侧落成首个 Monte Carlo acceptance-band test。
- `confirmed`: `_automation/test-artifacts/parity/e8_04_clustering_layer4_monte_carlo_contract.yaml`
  的 `coverage_large_g` contract 已被
  `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`
  直接消费，而不是手写阈值。
- `confirmed`: 当前 R helper 对该场景得到 `coverage_rate = 0.91`，
  位于 contract 要求的 `[0.90, 0.99]` band 内。

## 代码与测试

- 新增内部 helper：`lwdid-r/R/clustering_diagnostics.R`
  - `.parse_normal_distribution()`
  - `.draw_cluster_sizes()`
  - `.draw_cluster_treatment()`
  - `.simulate_clustering_monte_carlo_panel()`
  - `.run_clustering_monte_carlo_scenario()`
- 新增 Layer 4 regression：
  `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`
  - `E8-04 Layer 4 parity freezes standard coverage_large_g acceptance band`

## 验证证据

- RED:
  `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics-parity', reporter='summary')"`
  曾明确失败于缺失 `.run_clustering_monte_carlo_scenario()`
- GREEN:
  `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics', reporter='summary')"`
  当前通过
- GREEN:
  `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics-parity', reporter='summary')"`
  当前通过

## 剩余缺口

- Layer 4 仍剩 8 个 scenarios 待接入：
  `coverage_medium_g`、`coverage_small_g`、
  `wild_bootstrap_coverage_small_g`、`wild_bootstrap_relative_gap`、
  `size_under_null`、`power_under_alternative`、
  `unbalanced_treatment_coverage`、
  `heterogeneous_cluster_sizes_coverage`
- 其中 wild-bootstrap 两个场景尚未接入 R 侧 comparator；当前 helper 只覆盖
  standard cluster-robust path。
