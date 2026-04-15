# story-E8-04 Layer 3 smoking common-timing parity

## 时间

- `2026-03-24 12:42 CST`

## 本轮目标

- 为 `story-E8-04` 补齐 common-timing / built-in empirical 的 supplemental
  Layer 3 public-workflow parity。
- 直接复用 Python `smoking.csv` 与 R 包内置 `smoking` 数据，冻结
  `diagnose_clustering() -> recommend_clustering() -> check_clustering_consistency()`
  在 real-data common-timing 路径上的对齐结果。

## 新增工件

- `_automation/test-artifacts/parity/e8_04_clustering_layer3_smoking_common_timing_oracle.py`
- `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer3-smoking-common-timing.json`
- `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`

## Source anchors

- `lwdid-py_v0.2.3/data/smoking.csv`
- `lwdid-r/data/smoking.rda`
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_to_stata.py::common_timing_data`
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_to_stata_numerical.py::smoking_data`
- `_automation/test-artifacts/parity/20260324-theory-parity-e8-04-layer34-source-contract.md`

## 关键结论

1. 在 `smoking` common-timing real-data 路径上，Python 与当前 R 都把 `state`
   识别为：
   - `treatment_variation_level`
   - `recommended_var`
   - `check_clustering_consistency()` 的一致聚类层级
2. 该场景当前不需要 Wild Cluster Bootstrap：
   - `n_clusters = 39`
   - `n_treated_clusters = 1`
   - `n_control_clusters = 38`
   - `use_wild_bootstrap = FALSE`
3. `check_clustering_consistency(cluster_var = "state")` 当前继续保持：
   - `cluster_level = "same"`
   - `pct_inconsistent = 0`
   - `recommendation = "Clustering choice is appropriate."`
4. 这份证据属于 built-in empirical / common-timing supplemental Layer 3，
   不是 archived clustering tests 的 primary truth source；其作用是把
   `story-E8-04` 尚未冻结的 real-data common-timing 路径补成 executable oracle。

## TDD 轨迹

- RED：
  `test-clustering-diagnostics-parity.R` 先新增
  `E8-04 Layer 3 parity freezes smoking common-timing public workflow`，
  初次 rerun 失败，直接报缺
  `20260324-story-worker-e8-04-layer3-smoking-common-timing.json`。
- GREEN：
  新增 oracle generator 产出 JSON 后，重新执行
  `devtools::test(filter="clustering-diagnostics-parity")` 通过。
- VERIFY：
  再执行
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics", reporter="summary")'`
  当前通过，确认未破坏既有 clustering regression。

## 边界说明

- 本轮补的是 common-timing / built-in empirical supplemental Layer 3，
  不把它表述成论文直接给出的聚类阈值真值。
- 本轮未修改任何 tolerance，也未发现新的 Python-paper 冲突，因此
  `Docs/Python包bug列表.md` 保持不变。
- 结合已存在的 state-policy / few-cluster / hierarchical 三条 Layer 3 证据，
  `story-E8-04` 的直接剩余 gap 现只剩 Layer 4 Monte Carlo acceptance-band
  parity。
