# E8-04 Layer 4 wild-bootstrap scenarios

## 本轮结论

- `confirmed`: `wild_bootstrap_coverage_small_g` 与
  `wild_bootstrap_relative_gap` 已被 fresh rerun 固定为 executable Layer 4
  regression。
- `confirmed`: 两条场景都直接读取
  `_automation/test-artifacts/parity/e8_04_clustering_layer4_monte_carlo_contract.yaml`
  的 seed / simulation-count / WCB contract，而不是手写阈值。
- `confirmed`: Python default WCB contract 当前已被 R 侧 helper 复现为
  `requested_n_bootstrap = 199`、`actual_n_bootstrap = 1024`、
  `weight_type = rademacher`、`impose_null = TRUE`、
  `ci_method = percentile_t` 与 `full_enumeration = TRUE`。`story-E8-04`
  当前已无剩余 parity blocker。

## 数值证据

- `wild_bootstrap_coverage_small_g`:
  `coverage_rate = 0.96`，contract band `(0.85, +inf)`
- `wild_bootstrap_relative_gap`:
  `coverage_std_rate = 0.91`、`coverage_wild_rate = 0.96`、
  `coverage_rate_gap = 0.05`，contract band `[-0.10, +inf)`

## 验证命令

- GREEN:
  `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics-parity', reporter='summary')"`
- SUPPORT:
  `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics|clustering-diagnostics-parity', reporter='summary')"`

## 边界裁决

- `requested_n_bootstrap = 199` 只是 Python archived tests 的请求值，不是最终执行次数。
- 在 `G = 10` 且 `weight_type = rademacher` 时，Python default WCB contract 会自动转成
  `2^10 = 1024` 个权重组合的 full enumeration；R 当前结果已与该 effective contract 对齐。
- 该结论属于 Python-backed QA evidence，不应写成 `Docs/lw2026.md`
  的 paper theorem。

## 下一步建议

- 由 controller fresh rerun clustering suites 并复核本文件 / 同名 JSON。
- 若复核无误，应把 active story 从 `story-E8-04` 切换到 `story-E8-05`。
