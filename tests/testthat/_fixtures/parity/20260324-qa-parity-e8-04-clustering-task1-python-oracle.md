# story-E8-04 Task 1 Python oracle parity

## 时间

- `2026-03-24 06:58 CST`

## 本轮目标

- 为 `story-E8-04` 已落地的 `Task E8-04.1` 建立 Python-backed shared fixture /
  oracle / testthat 回归，避免当前聚类诊断 helper 继续只靠手写常数断言。

## 新增工件

- `_automation/test-artifacts/parity/e8_04_clustering_task1_oracle.py`
- `_automation/test-artifacts/parity/e8_04_clustering_hierarchical_fixture.csv`
- `_automation/test-artifacts/parity/e8_04_clustering_small_cluster_fixture.csv`
- `_automation/test-artifacts/parity/e8_04_clustering_within_cluster_variation_fixture.csv`
- `_automation/test-artifacts/parity/e8_04_clustering_never_treated_fixture.csv`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-clustering-task1-python-oracle.json`
- `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`

## Source anchors

- Python archived tests:
  `lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py`
- Python implementation:
  `lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py`
- R implementation:
  `lwdid-r/R/clustering_diagnostics.R`

## 关键发现

1. `parity-qa-status.yaml` 中“`clustering_diagnostics.R` / `test-clustering-diagnostics.R`
   仍不存在”的 blocker 已过时；这两个文件已在 `2026-03-24 06:43 CST` 被
   `story-worker` 落地。
2. 新增的 Python oracle 现已覆盖三类 helper-level synthetic cases：
   - archived hierarchical fixture：锁定 `higher/same` 层级判定与
     `state/county` 的 `ClusterVarStats`
   - archived within-cluster-variation fixture：锁定
     `treatment_varies_within`
   - QA never-treated fixture：锁定 `gvar = NA/0/Inf` 的 treated-mask 语义
3. fresh RED/GREEN parity 回归额外抓出一个真实 drift：
   Python `groupby(...).nunique()` 默认忽略 `NA`，而 R 原实现用
   `data.table::uniqueN()` 把 `NA` 计为额外 unique value，导致
   `never-treated (NA/0)` 聚类被误判为 `within-cluster variation`。
   本轮已将 R 侧对齐为 `uniqueN(..., na.rm = TRUE)`，并同步修正
   `test-clustering-diagnostics.R` 的旧常数。

## 验证

- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics", reporter="summary")'`
  当前通过。
- 新增 `test-clustering-diagnostics-parity.R` 已实际读取
  `20260324-qa-parity-e8-04-clustering-task1-python-oracle.json`，并对比 R / Python
  的 helper-level exact / numeric 结果。

## 边界更新

- `Task E8-04.1` 现在具备 executable Layer 1 / Layer 2 parity 证据，但仅限
  helper-level synthetic fixtures。
- 当前 story blocker 不再是“缺文件”，而是 `Task E8-04.2` / `Task E8-04.3`
  尚未实现，full-module exact / numeric parity 仍不能 closure-ready。
- 本轮未变更任何 tolerance，也未新增 Python bug ledger 条目。
