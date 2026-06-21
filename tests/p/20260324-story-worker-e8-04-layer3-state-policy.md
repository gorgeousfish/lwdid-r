# story-E8-04 Layer 3 state-policy parity

## 时间

- `2026-03-24 08:47 CST`

## 本轮目标

- 为 `story-E8-04` 启动首个 Layer 3 real-data-like parity regression。
- 直接复用 archived Python clustering scenario，把 state-policy / county-observation
  面板上的 public workflow 冻结成 R 侧可执行 oracle。

## 新增工件

- `_automation/test-artifacts/parity/e8_04_clustering_layer3_state_policy_oracle.py`
- `_automation/test-artifacts/parity/e8_04_clustering_layer3_state_policy_fixture.csv`
- `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer3-state-policy.json`
- `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`

## Source anchors

- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_simulated.py`
  - `TestStateLevelPolicy::test_diagnose_recommends_state`
  - `TestStateLevelPolicy::test_recommend_no_wild_bootstrap`
  - `TestStateLevelPolicy::test_consistency_check_passes`
  - `TestStateLevelPolicy::test_state_is_higher_level`
- `_automation/test-artifacts/parity/20260324-theory-parity-e8-04-layer34-source-contract.md`

## 关键结论

1. 在 50-state / 1000-county / 10-period 的 state-policy fixture 上，
   Python 与当前 R 都把 `state` 识别为：
   - `treatment_variation_level`
   - `recommended_cluster_var`
   - `check_clustering_consistency()` 的一致聚类层级
2. 该场景下 `recommend_clustering()` 当前不需要 Wild Cluster Bootstrap：
   - `n_clusters = 50`
   - `use_wild_bootstrap = FALSE`
   - 仍保留 `county` 作为 alternative，并把其理由冻结为
     `reliability score: 1.00`
3. 这份证据说明 `story-E8-04` 的 Layer 3 已经从 “未启动” 进入
   “已起步但未收口” 状态。剩余 gap 不再是第一条 real-data-like regression，
   而是 few-cluster / hierarchical supplemental Layer 3 与 Layer 4 Monte Carlo。

## TDD 轨迹

- RED：
  `test-clustering-diagnostics-parity.R` 先新增
  `E8-04 Layer 3 parity freezes state-policy public workflow`，
  初次 rerun 失败，直接报缺
  `20260324-story-worker-e8-04-layer3-state-policy.json`。
- GREEN：
  新增 oracle generator 产出 fixture + JSON 后，重新执行
  `devtools::test(filter="clustering-diagnostics-parity")` 通过。
- VERIFY：
  再执行
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics", reporter="summary")'`
  当前通过，确认未破坏现有 clustering suite。

## 边界说明

- 本轮只启动了 Layer 3 的 state-policy regression，不声称 `Task E8-04.5` 已完成。
- 本轮未修改任何 tolerance，也未发现新的 Python-paper 冲突，因此
  `Docs/Python包bug列表.md` 保持不变。
