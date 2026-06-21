# E8-04 Layer 4 standard cluster-robust scenarios

## 本轮结论

- `confirmed`: `coverage_medium_g`、`coverage_small_g`、`size_under_null`、
  `power_under_alternative`、`unbalanced_treatment_coverage` 与
  `heterogeneous_cluster_sizes_coverage` 已被 fresh rerun 固定为 executable
  Layer 4 regression。
- `confirmed`: 上述 6 个场景全部直接读取
  `_automation/test-artifacts/parity/e8_04_clustering_layer4_monte_carlo_contract.yaml`
  的 seed / simulation-count / acceptance-band 设定，而不是手写阈值。
- `confirmed`: 结合既有
  `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer4-coverage-large-g.json/.md`，
  Layer 4 的 standard cluster-robust path 已全部冻结；当前剩余 blocker 只剩
  `wild_bootstrap_coverage_small_g` 与 `wild_bootstrap_relative_gap`。

## 数值证据

- `coverage_medium_g`: `coverage_rate = 0.92`，contract band `[0.88, 0.99]`
- `coverage_small_g`: `coverage_rate = 0.93`，contract band `(0.80, +inf)`
- `size_under_null`: `rejection_rate = 0.065`，contract band `[0.02, 0.10]`
- `power_under_alternative`: `power = 0.72`，contract band `(0.50, +inf)`
- `unbalanced_treatment_coverage`: `coverage_rate = 0.91`，contract band `(0.85, +inf)`
- `heterogeneous_cluster_sizes_coverage`: `coverage_rate = 0.92`，contract band `(0.85, +inf)`

## 验证命令

- RED:
  `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics-parity', reporter='summary')"`
  先前失败于缺失
  `20260324-qa-parity-e8-04-layer4-standard-scenarios.json`
- GREEN:
  `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics', reporter='summary')"`
- GREEN:
  `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics-parity', reporter='summary')"`

## 剩余缺口

- `wild_bootstrap_coverage_small_g`
- `wild_bootstrap_relative_gap`
- 当前 R helper 仍只覆盖 standard cluster-robust path；若下一轮继续推进，
  应优先决定是最小实现 wild cluster bootstrap CI，还是用已有外部实现做 QA
  bridge，而不是回退去重写 standard scenarios。
