# E8-04 Layer 3 / Layer 4 source contract

## 本轮结论

`story-E8-04` 当前的理论帮工重点已经从 public API 缺口，转移为
Layer 3 / Layer 4 parity 的合法来源界定。本轮 source audit 结论如下：

1. `Docs/lw2026.md:641-643` 只提供聚类推断的论文真值边界：
   当处理在高于观测单位的层级分配，且 treated / control clusters 数量足够时，
   可以在更高层级聚类；该论证建立在 `lwdid` 变换后的横截面回归表示上。
2. Python archived clustering tests 才是 `story-E8-04` 下一层 parity 的直接可执行来源。
   它们提供了 real-data-like hierarchical scenarios、few-cluster WCB 场景，
   以及 Monte Carlo coverage / size / power acceptance bands。
3. 这些 Layer 3 / 4 阈值属于 Python QA contract 与项目设计扩展，
   不是论文直接给出的定理阈值；因此后续 R parity 不应把它们表述成
   “paper theorem”，而应表述为 “Python-backed QA target”。
4. 本轮未发现新的 Python-paper 冲突，因此不更新 `Docs/Python包bug列表.md`。

## Source anchors

### 1. 论文真值

- `Docs/lw2026.md:641-643` 明确指出：`lwdid` 经过时间序列变换后，
  推断基于横截面回归；若政策只在州级变化而观测单位是县，则可在州级聚类，
  前提是 treated / control states 数量足够。
- 论文没有给出 `G >= 20`、`coverage >= 0.90`、`wild bootstrap > 0.85`
  这类数值阈值；这些都不是 paper-only truth。

### 2. Python Layer 3 场景来源

- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_simulated.py:34-95`
  生成 state-policy / county-observation 场景：
  `n_states = 50`、`counties_per_state = 20`，并要求推荐 `state`、
  `n_clusters = 50`、`use_wild_bootstrap = FALSE`、一致性检查通过。
- 同文件 `:200-262` 生成 few-cluster 场景：
  `n_regions = 5`，要求 `use_wild_bootstrap = TRUE`，
  且 `wild_cluster_bootstrap(..., n_bootstrap = 99)` 可以返回合法 `pvalue`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_empirical.py:300-351`
  再次固定 hierarchical panel / small-cluster 的 public contract：
  `state` 为推荐层级，small-cluster 场景下需要 WCB，
  且存在 alternatives 列表。

### 3. Python Layer 4 场景来源

- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:29-97`
  固定 cluster-robust CI coverage 的 acceptance bands：
  - `G = 50`: `0.90 <= coverage_rate <= 0.99`
  - `G = 20`: `0.88 <= coverage_rate <= 0.99`
- 同文件 `:100-134` 固定 small-cluster baseline：
  - `G = 10`: `coverage_rate > 0.80`
- 同文件 `:141-230` 固定 WCB contract：
  - `G = 10` 时 `wild bootstrap coverage > 0.85`
  - 且 `coverage_wild_rate >= coverage_std_rate - 0.10`
- 同文件 `:248-352` 进一步固定 Monte Carlo 的 size / power / robustness bands：
  - `G = 30`、`tau = 0` 时拒绝率 `0.02 <= rejection_rate <= 0.10`
  - `G = 30`、`tau = 2` 时 `power > 0.50`
  - 不平衡处理分配时 `coverage_rate > 0.85`
  - 异质 cluster sizes 时 `coverage_rate > 0.85`

## 对 E8-04 的直接合同

### Layer 3 应该测什么

- 应优先复用 Python archived tests 中已经成型的 clustering scenarios：
  state-policy / county、few clusters、hierarchical panel、multi-level structure。
- 可执行 parity 的核心输出应是：
  - `recommended_cluster_var`
  - `n_clusters`
  - `use_wild_bootstrap`
  - `pct_inconsistent`
  - `level_relative_to_unit`
  - `treatment_variation_level`
- 若后续复用 `smoking` / `castle`，只能把它们当作 supplemental scenario，
  不能把它们写成 clustering 模块唯一的 source-backed Layer 3 真值；
  因为 Python clustering archived tests 当前并不是围绕这两个内置数据集定义的。

### Layer 4 应该怎么验

- Monte Carlo parity 不应追求 exact-number oracle；正确合同是 band-based regression。
- 推荐直接沿用 Python archived tests 的 acceptance bands，而不是重新拍脑袋设阈值。
- 由于随机模拟存在 Monte Carlo noise，R 端若做 parity，应该固定：
  - seed
  - `n_simulations`
  - `obs_per_cluster`
  - `n_bootstrap`
  然后对 coverage / rejection / power 使用区间断言，而不是 exact equality。

## 非 claim 边界

- `G >= 20` 与 WCB 建议属于 Python / 项目设计中的 operational rule，
  受 Cameron, Gelbach & Miller (2008) 及 Abadie et al. (2023) 的启发，
  但不是 `Docs/lw2026.md` 的逐字定理阈值。
- Layer 4 coverage / size / power 阈值属于 QA acceptance bands，
  不应在发布文档中表述成论文数值结论。
- 现阶段没有证据支持把 `story-E8-04` 写成 Layer 3 / 4 已完成；
  本轮只锁定了这些层的 comparator contract。

## 给下一个角色的具体建议

- `qa-parity`：
  先把 `test_clustering_simulated.py` 的 state-policy / few-cluster 场景抽成
  R fixture + JSON oracle，再决定是否扩展到 multi-level / empirical panel。
- `story-worker`：
  若需要补 `Task E8-04.4/5` 测试，不要把 Layer 4 写成 exact-number test；
  应直接采用本 note 固定的 acceptance bands。
- `controller`：
  在 `story-E8-04` 未把 Layer 3 / 4 regression 落地前，不应把
  `tasks.md` 的 closure-ready 目标语句误读为“当前已完成”。
