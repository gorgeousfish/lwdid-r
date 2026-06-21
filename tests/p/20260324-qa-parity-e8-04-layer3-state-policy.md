## E8-04 Layer 3 state-policy public workflow

- 时间：2026-03-24 08:51 CST
- 角色：qa-parity
- phase：implementation
- active_story：story-E8-04

### 本轮执行

先重放 `_automation/test-artifacts/parity/e8_04_clustering_layer3_state_policy_oracle.py`，
确认 `20260324-story-worker-e8-04-layer3-state-policy.json` 仍由当前
`lwdid-py_v0.2.3` 生成；随后 fresh rerun
`Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics-parity", reporter="summary")'`
与
`Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics", reporter="summary")'`。

### 关键证据

- `e8_04_clustering_layer3_state_policy_fixture.csv` 当前固定为 50 个 state、
  1000 个 county、10 个时期的 state-level treatment policy panel。
- `20260324-story-worker-e8-04-layer3-state-policy.json` 当前继续给出：
  `recommended_cluster_var = "state"`、
  `treatment_variation_level = "state"`、
  `python_consistency$is_consistent = TRUE`、
  `pct_clusters_with_variation = 0`。
- fresh rerun `clustering-diagnostics-parity` 当前通过，说明 R 侧
  `diagnose_clustering()`、`recommend_clustering()` 与
  `check_clustering_consistency()` 已能直接执行该 Layer 3 state-policy public
  workflow，而不再停留在一次性 JSON。
- fresh rerun `clustering-diagnostics` 当前也通过，说明该 Layer 3 contract
  没有把已有的 Task 1/2/3 与 Task 4 regression 打回红色。

### 结论

`story-E8-04` 的 Layer 3 已不再是“尚未启动”。当前至少已有一个
Python-backed state-policy public-workflow scenario 具备 executable parity。
剩余缺口转为：

- few-cluster / built-in empirical clustering scenarios 尚未冻结为同级 oracle
- Layer 4 Monte Carlo acceptance-band parity 尚未落地

### Python bug 审计

本轮未发现新的 Python-paper drift，不更新 `Docs/Python包bug列表.md`。
