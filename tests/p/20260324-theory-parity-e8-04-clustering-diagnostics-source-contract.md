# story-E8-04 clustering diagnostics source contract

## 时间

- `2026-03-24 06:29 CST`

## Target

- 为 `story-E8-04` 冻结聚类诊断模块的 source-backed 真值合同，避免
  `story-worker` 在 `R/clustering_diagnostics.R` 尚不存在时误抄 Python 公共接口
  的排序细节或把经验阈值写成论文定理。
- 明确哪些规则来自论文与计量原理，哪些来自 Python 参考实现，哪些是当前 R
  spec 的有意设计差异。

## Highest-Truth Anchors

- `Docs/lw2026.md:641-643`：Lee-Wooldridge 2026 明确把变换后的估计写成
  cross-sectional regression，并指出当政策在高于观测单位的层级变化时，
  只要 treated/control clusters 数量足够，就应在该更高层级聚类。
- `Docs/lw2026.md:25-31`, `Docs/lw2026.md:651-657`：论文把小样本推断主线放在
  cross-sectional exact / HC3 inference 上；`G < 20` 时建议 WCB 属于
  CGM (2008) / Abadie et al. (2023) 的实践规则，而不是 lw2026 的显式定理。
- `.kiro/specs/story-E8-04/requirements.md:48-63`：E8-04 已把 score、排序、
  never-treated 识别和 consistency 容差固化为实现契约。

## Source-Backed Contract

### 1. Paper-backed rule

- 真正需要自动化的是“在哪个更高层级聚类”，而不是改写 ATT 公式本身。
- 因为 `lwdid` 先把面板数据变成截面回归，聚类标准误选择应围绕“处理在哪个层级变异”
  来做；若单位是 county、政策在 state 变异，则应优先考虑 state 聚类。
- 论文没有给出 `0.5/0.3/0.2` 的 reliability score，也没有给出 `G < 20`
  的 WCB 阈值；这些属于实现层 heuristic / literature-backed rule，
  不得误表述为 lw2026 的公式。

### 2. Python-backed arithmetic and comparator expectations

- `ClusterVarStats.reliability_score`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:185-214`) 定义为
  `0.5 * min(G / 50, 1) + 0.3 * balance_score + 0.2 * cv_score`。
- `balance_score`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:203-208`) 使用
  `min(n_treated_clusters, n_control_clusters) / (n_clusters / 2)`，再截断到 `1.0`。
- `cv_score`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:210-214`) 使用
  `max(0, 1 - cluster_size_cv / 2)`。
- `cluster_sizes`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:598-600`) 基于
  `data.groupby(cluster_var).size()`，因此真值对象是“观测数”，不是去重后的 unit 数。
- gvar 模式 treated mask
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:602-606`) 将
  `0` / `Inf` / `NA` 视为 never treated；R 侧必须复制这个 treated-mask 逻辑。
- `treatment_varies_within_cluster`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:631-640`) 用
  `groupby(cluster_var)[gvar_or_d].nunique()` 判断，关键对象是“唯一值个数”，
  不是 treated/control 是否共存；且 pandas 默认 `dropna=True`，所以
  `gvar = {NA, 0}` 的 never-treated cluster 不算 within-cluster variation，
  而 `gvar = {Inf, 4}` 仍算 variation。
- `detect_treatment_variation_level`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:697-710`) 先按
  `nunique()` 升序遍历候选变量，再返回第一个“组内 treatment unique count 恒为 1”
  的层级；若所有层级都变，则回落到 `ivar`。
- `_generate_clustering_recommendation`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:734-782`) 的关键契约是：
  - 先过滤 `is_valid_cluster`
  - treatment-level 若 `G >= 20` 立即推荐
  - treatment-level 若 `G < 20` 只发 warning，不提前 return
  - 回退排序使用双键 `(n_clusters >= 20, reliability_score)`
- `check_clustering_consistency`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:1272-1289`) 的一致性判定是：
  `pct_with_variation < 5` 且 `cluster_level in {'same', 'higher'}`。

### 3. Intentional R design deltas

- `recommend_clustering_level()` 的 Python 公共接口
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:1144-1164`) 只按
  `reliability_score` 单键排序。
- 但 `.kiro/specs/story-E8-04/design.md:259` 明确规定：
  R 的 `recommend_clustering()` 要复用 `diagnose_clustering()` 的双层排序
  `(n_clusters >= 20, reliability_score)`。
- 该差异是 `design-delta`，不是 `Python bug`。后续 parity / review 不应要求
  R public recommender 回退到 Python 的单键排序。
- 同理，R public 接口在“无有效聚类选项”时返回 `recommended_var = NULL`
  是 spec 明确的 R 惯例差异，不应视为 drift。

### 4. Consistency-check scope note

- Python `check_clustering_consistency()`
  (`lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:1278-1280`) 调用
  `_detect_treatment_variation_level(data, ivar, [cluster_var], ...)`，
  只把当前 `cluster_var` 作为候选集合传入。
- 因此该函数返回的 `treatment_variation_level` 是“相对于当前被检查变量的局部判断”，
  不是对所有候选层级做一次全局搜索。
- 若 R 侧想保留 Python 语义，应维持这个局部范围；若未来要改成全局搜索，
  必须作为显式设计变更记录，不能默默改写。

## Bug Judgment

- `confirmed-no-new-bug`: 本轮未发现新的 Python-paper 冲突，不更新
  `Docs/Python包bug列表.md`。
- `confirmed`: `recommend_clustering_level()` 的单键排序属于 Python 公共接口行为，
  但已被 E8-04 spec 显式覆盖为 R 设计差异；不得误记为 bug。

## Implementation Implications

1. `story-worker` 应先把 `REQ-11` 到 `REQ-25` 直接翻译为 helper-level tests，
   再实现 `R/clustering_diagnostics.R`。
2. `TC-8.4.1` / `TC-8.4.13` 的数值真值应按 score component 分别锁定，
   不要把 `WCB` warning 文案与 reliability arithmetic 混在一个断言里。
3. `TC-8.4.8` / `TC-8.4.9` 必须明确区分：
   - internal recommender 的双键排序；
   - public `recommend_clustering()` 的有意 R 设计差异；
   - Python public recommender 的单键排序。
4. 当前仓内已落地 `lwdid-r/R/clustering_diagnostics.R`、
   `lwdid-r/tests/testthat/test-clustering-diagnostics.R` 与
   `test-clustering-diagnostics-parity.R`；后续重点应转向
   `Task E8-04.2 / 3` 的 helper/public API，而不是继续把 “Task 1 缺文件”
   当作 blocker。
