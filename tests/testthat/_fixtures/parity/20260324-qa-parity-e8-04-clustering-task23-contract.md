# story-E8-04 Task 3 public-contract oracle

## 时间

- `2026-03-24 07:50 CST`

## 本轮目标

- 在不越过 `qa-parity` 角色边界的前提下，为 `story-E8-04` 的
  `Task E8-04.3` 冻结可执行的 public-contract oracle，重点覆盖：
  - small-cluster 场景下的 public diagnosis / recommendation 语义；
  - `check_clustering_consistency()` 的 local-scope +
    `nunique(dropna=True)` 语义。

## 新增工件

- `_automation/test-artifacts/parity/e8_04_clustering_task23_contract_oracle.py`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-clustering-task23-contract.json`
- `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`

## Source anchors

- Python archived tests:
  `lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py`
- Python implementation:
  `lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py`
- theory-parity contract:
  `_automation/test-artifacts/parity/20260324-theory-parity-e8-04-consistency-na-scope-contract.md`

## 关键发现

1. 复用现有 shared fixtures 后，small-cluster 场景的 Python public outputs 已可被
   单独冻结为 oracle：
   - `diagnose_clustering()` 仍推荐 `state`；
   - `recommend_clustering_level()` 仍要求 `use_wild_bootstrap = TRUE`；
   - 当前 public recommendation 的关键锚点是 `G = 5`、`treated = 2`、
     `control = 3`、`confidence = 0.49`。
2. `check_clustering_consistency()` 的 local-scope 合同现已被可执行 JSON 固化：
   - hierarchical fixture 上，`cluster_var = state` 时
     `pct_clusters_with_variation = 0` 且 `treatment_variation_level = state`；
   - never-treated fixture 上，`cluster_var = cluster` 时
     `cluster_unique_counts = {A:1, B:1, C:2}`，
     `pct_clusters_with_variation = 33.333...%`，
     `treatment_variation_level = id`。
3. 这组证据再次证明 `dropna=True` 是必须保留的 API 语义：
   `{NA, 0}` 只能算 single-valued never-treated cluster，不能被误放大成
   `66.666...%` inconsistency。

## 验证

- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics-parity", reporter="summary")'`
  当前通过。
- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics", reporter="summary")'`
  当前通过。

## 边界更新

- `Task E8-04.2` 的 `.detect_treatment_variation_level()` 与
  `Task E8-04.3` 的 consistency public contract 现在都已有 executable oracle；
  parity blocker 不应再写成“helper/public contract 未冻结”。
- 当前剩余 story-level gap 已进一步收敛为 R public API
  `recommend_clustering()` / `check_clustering_consistency()` 及对应 print/S3
  输出尚未落地，而不是 fixture / oracle 缺失。
- 本轮未变更任何 tolerance，也未新增 Python bug ledger 条目。
