## E8-04 Layer 3 hierarchical public workflow

- 时间：2026-03-24 12:33 CST
- 角色：qa-parity
- phase：implementation
- active_story：story-E8-04

### 本轮执行

先在 `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R` 新增
hierarchical Layer 3 RED 用例，并 fresh rerun
`Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics-parity", reporter="summary")'`；
失败原因精确定位为缺失
`_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer3-hierarchical.json`。

随后新增
`_automation/test-artifacts/parity/e8_04_clustering_layer3_hierarchical_oracle.py`，
按 archived Python `test_clustering_empirical.py` 的
`hierarchical_panel_data` / `TestIntegrationEmpirical::test_full_workflow`
工作流重建 shared fixture
`e8_04_clustering_layer3_hierarchical_fixture.csv` 与 JSON oracle
`20260324-qa-parity-e8-04-layer3-hierarchical.json`，再 fresh rerun
`devtools::test(filter="clustering-diagnostics|clustering-diagnostics-parity")`
确认新旧 Layer 3 regression 同时通过。

### 关键证据

- `e8_04_clustering_layer3_hierarchical_fixture.csv` 当前固定为 4 个 region、
  20 个 state、200 个 industry、1000 个 idcode、10 个时期的层级 panel；
  `first_treat` 在 `state < 10` 时为 `6`，否则为 `0`。
- `20260324-qa-parity-e8-04-layer3-hierarchical.json` 当前给出：
  `python_diagnose_clustering$recommended_cluster_var = "state"`、
  `python_diagnose_clustering$treatment_variation_level = "region"`、
  `python_recommend_clustering$recommended_var = "state"`、
  `python_recommend_clustering$use_wild_bootstrap = FALSE`、
  `python_consistency$is_consistent = TRUE`、
  `python_consistency$pct_clusters_with_variation = 0`。
- fresh rerun
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics-parity", reporter="summary")'`
  当前通过，说明 hierarchical `diagnose -> recommend -> consistency`
  public workflow 已具备 executable Layer 3 parity。
- fresh rerun
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics|clustering-diagnostics-parity", reporter="summary")'`
  当前通过，说明新增 hierarchical oracle 没有把已有
  state-policy / few-cluster / Task 4 / Task 5 回归打回红色。

### 结论

`story-E8-04` 当前已具备三个可执行 Layer 3 public-workflow scenario：

- state-policy
- few-cluster
- hierarchical panel

剩余缺口进一步收敛为：

- built-in empirical / common-timing supplemental Layer 3 场景尚未冻结
- Layer 4 Monte Carlo acceptance-band parity 尚未落地

### Python bug 审计

本轮未发现新的 Python-paper drift，不更新 `Docs/Python包bug列表.md`。
