# E8-04 clustering public API state convergence

- 时间：`2026-03-24 08:04 CST`
- 角色：`correct-course`
- active story：`story-E8-04`

## 本轮结论

`story-E8-04` 先前的真实 drift 是：spec / parity 状态已经把 `Task E8-04.3`
收紧到 public API contract，但 R 代码事实仍缺
`recommend_clustering()`、`check_clustering_consistency()` 与 clustering
print/S3/NAMESPACE 合同。该 drift 现已解除。

## 本轮证据

1. `lwdid-r/tests/testthat/test-clustering-diagnostics.R` 已新增 public API RED→GREEN 回归：
   - `recommend_clustering()` 返回 `lwdid_clustering_recommendation`
   - `check_clustering_consistency()` 返回 9 字段 contract
   - generic `print()` 现触发 clustering diagnosis / recommendation 的 S3 输出
   - `getNamespaceExports("lwdid")` 现包含 `recommend_clustering` 与 `check_clustering_consistency`
2. `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R` 已把
   `20260324-qa-parity-e8-04-clustering-task23-contract.json` 接到 public API：
   - `recommend_clustering()` 现对齐 small-cluster `recommended_var = "state"`、
     `confidence = 0.49`、`use_wild_bootstrap = TRUE` 与 Python reasons / warnings
   - `check_clustering_consistency()` 现对齐 hierarchical `0%` 与 never-treated
     fixture `33.333...%` 的 inconsistency share、details 与 recommendation
3. `lwdid-r/R/clustering_diagnostics.R` 已落地：
   - `.generate_recommendation_reasons()` / `.get_alternative_reason()`
   - `recommend_clustering()`
   - `check_clustering_consistency()`
   - `print.lwdid_clustering_diagnosis()`
   - `print.lwdid_clustering_recommendation()`
4. `Rscript -e 'devtools::document(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r")'`
   已把 `lwdid-r/NAMESPACE` 的 export / S3 注册同步到代码事实。
5. fresh rerun：

```r
Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics", reporter="summary")'
```

当前通过。

## 状态收敛

- 已解除的 drift：
  - `Task E8-04.2 / 3` code-vs-spec drift
  - clustering public API 缺失导致的 state-machine drift
  - clustering print/S3/NAMESPACE contract 缺口
- 仍未解除的 gap：
  - Task `4/5` 的扩展测试矩阵
  - Layer 3 / Layer 4 的 real-data / Monte Carlo parity
  - `story-E8-05` 仍保持下一顺位，不得抢占 active slot
