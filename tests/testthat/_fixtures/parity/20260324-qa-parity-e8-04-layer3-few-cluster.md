## E8-04 Layer 3 few-cluster public workflow

- 时间：2026-03-24 10:07 CST
- 角色：qa-parity
- phase：implementation
- active_story：story-E8-04

### 本轮执行

先在 `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R` 新增
few-cluster Layer 3 RED 用例，并 fresh rerun
`Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics-parity", reporter="summary")'`；
失败原因精确定位为缺失
`_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer3-few-cluster.json`。

随后新增
`_automation/test-artifacts/parity/e8_04_clustering_layer3_few_cluster_oracle.py`，
按 archived Python `TestFewClusters` 场景重建 shared fixture
`e8_04_clustering_layer3_few_cluster_fixture.csv` 与 JSON oracle
`20260324-qa-parity-e8-04-layer3-few-cluster.json`。同一轮 green rerun 还暴露出
R 端 `diagnose_clustering()` 在 `<20 clusters` 路径上的 warning string
仍未与 Python 对齐；现已在 `lwdid-r/R/clustering_diagnostics.R` 与对应
`test-clustering-diagnostics*.R` 中收口为 Python-backed contract。

### 关键证据

- `e8_04_clustering_layer3_few_cluster_fixture.csv` 当前固定为 5 个 region、
  1000 个 unit、8 个时期的 few-cluster panel，处理在 region 层级分配，
  `first_treat ∈ {5, 0}`。
- `20260324-qa-parity-e8-04-layer3-few-cluster.json` 当前给出：
  `recommended_cluster_var = "region"`、
  `treatment_variation_level = "region"`、
  `python_recommend_clustering$use_wild_bootstrap = TRUE`、
  `python_consistency$is_consistent = TRUE`、
  `pct_clusters_with_variation = 0`。
- fresh rerun
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics-parity", reporter="summary")'`
  当前通过，说明 few-cluster `diagnose -> recommend -> consistency`
  public workflow 已具备 executable Layer 3 parity。
- fresh rerun
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics", reporter="summary")'`
  当前通过，说明 warning contract 对齐没有把现有 Task `1/2/3/4/5`
  regression 打回红色。

### 结论

`story-E8-04` 当前已不再只拥有 state-policy 一个 Layer 3 scenario；
few-cluster public workflow 现也具备 Python-backed executable evidence。
剩余缺口进一步收敛为：

- hierarchical / built-in empirical supplemental Layer 3 场景尚未冻结
- Layer 4 Monte Carlo acceptance-band parity 尚未落地
- Task `5.2` 的边界条件矩阵仍未补齐

### Python bug 审计

本轮未发现新的 Python-paper drift，不更新 `Docs/Python包bug列表.md`。
本轮发现的是 R 侧 warning string drift，而非 Python bug。
