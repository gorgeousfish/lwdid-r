# E8-04 Layer 4 standard path state convergence

## 本轮结论

- `confirmed`: `story-E8-04` 当前已不再是“Layer 4 Monte Carlo 仅有 contract、尚无
  standard path regression”的状态。
- `confirmed`: 结合既有
  `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer4-coverage-large-g.json/.md`
  与新增
  `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer4-standard-scenarios.json/.md`，
  Layer 4 的 standard cluster-robust 7 个 scenario 已全部具备
  executable regression。
- `confirmed`: 当前 canonical state 若仍把 `coverage_large_g`、
  `coverage_medium_g`、`coverage_small_g`、`size_under_null`、
  `power_under_alternative`、`unbalanced_treatment_coverage` 或
  `heterogeneous_cluster_sizes_coverage` 写为 pending，即构成 state drift。

## fresh verification

- `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics', reporter='summary')"`
- `Rscript -e "devtools::test(pkg='/Users/cxy/Desktop/lwdid_r/lwdid-r', filter='clustering-diagnostics-parity', reporter='summary')"`

以上两条命令在 `2026-03-24 14:00 CST` fresh rerun 下当前通过。

## 控制面裁决

- `story-E8-04` 仍保持 active story。
- Layer 4 当前剩余 blocker 只保留：
  - `wild_bootstrap_coverage_small_g`
  - `wild_bootstrap_relative_gap`
- `E8-05` 继续保持下一顺位，不得因旧的 Layer 4 宽泛文案而抢占 active slot。

## 容差与 bug ledger

- 本轮未放宽任何容差。
- 本轮未发现新的 Python-paper drift，因此 `Docs/Python包bug列表.md` 保持不变。
