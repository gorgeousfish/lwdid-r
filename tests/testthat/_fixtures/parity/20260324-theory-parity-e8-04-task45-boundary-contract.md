# story-E8-04 Task 4 / 5 boundary contract

## 时间

- `2026-03-24 09:36 CST`

## Target

- 为 `story-E8-04` 当前仍未收口的 `Task E8-04.4` / `Task E8-04.5`
  冻结 source-backed 边界合同，重点回答三个问题：
  1. 哪些断言属于 Python / 数学恒等式，可直接写成 exact regression；
  2. 哪些断言只是 R 设计差异，不应误记为 Python bug；
  3. `tasks.md` 草案里出现的 `G=0` 边界，当前是否已有可追随的 parity 真值。

## Highest-Truth Anchors

- 论文真值：
  `Docs/lw2026.md:641-643` 只给出“当处理在高于观测单位的层级分配时，可在该更高层级聚类”的横截面回归推断原则。
- Python helper / public API：
  `lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py`
  中的 `ClusterVarStats.reliability_score`、`is_valid_cluster`、
  `_generate_clustering_recommendation()`、
  `recommend_clustering_level()`、
  `check_clustering_consistency()`。
- 当前 spec：
  `.kiro/specs/story-E8-04/requirements.md` 的 `REQ-11` 到 `REQ-26`、
  `AC-02/04/05/08/09/10/12/13/15`；
  `.kiro/specs/story-E8-04/design.md` 中对 internal recommender 双键排序与
  public `recommend_clustering()` 设计差异的说明。

## Mainline Contract

### 1. `Task 4.1` 与 `Task 5.1` 的 exact arithmetic contract

- `reliability_score` 的唯一 source-backed 公式仍是：

  $$
  Score = 0.5 \times \min(G/50, 1) + 0.3 \times BalanceScore + 0.2 \times CVScore.
  $$

- 其中
  - `BalanceScore = min(min(G_t, G_c)/(G/2), 1)`，仅在 `G > 0` 时按该式计算；
    `G = 0` 时 Python / spec 都没有给出 paper-backed 数学意义。
  - `CVScore = max(0, 1 - CV/2)`。
- 因为这是纯代数合同，`Task 4.1` 与 `Task 5.1` 允许直接使用 exact numeric
  regression。以下值可作为 source-backed regression anchors：
  - `G=40, G_t=20, G_c=20, CV=0.5 -> 0.85`
  - `G=10, G_t=8, G_c=2, CV=1.5 -> 0.27`
  - `G=100, G_t=50, G_c=50, CV=0.1 -> 0.99`
  - `G=1, G_t=0, G_c=1, CV=0 -> 0.21`
  - `G=50, G_t=50, G_c=0, CV=0 -> 0.70`
  - `G=20, G_t=10, G_c=10, CV=3.0 -> 0.50`
- 这些值属于公式直接推导，不需要额外容差放宽；`1e-10` 级断言是合理的。

### 2. `Task 5.4` 的 BalanceScore 数学等价性

- 令 `p = G_t / G`，当 `G > 0` 时有

  $$
  \min\left(\frac{\min(G_t, G_c)}{G/2}, 1\right)
  = 2 \min(p, 1-p)
  = 1 - |2p - 1|
  = \max(0, 1 - 2|p - 0.5|).
  $$

- 因为 `p in [0, 1]`，上式中 `1 - |2p - 1|` 本身已非负；
  因此 `max(0, ...)` 只是把定义写成对边界更稳健的形式，不改变 `G > 0`
  的数学值。
- 这意味着 `Task 5.4` 应被视为恒等式验证，而不是经验近似测试。

### 3. `Task 4.3` / `Task 4.4` 的 strict inequality contract

- WCB 推荐阈值来自 Python public `recommend_clustering_level()` 的
  `best_stats.n_clusters < min_clusters`，以及 internal
  `_generate_clustering_recommendation()` 的 small-cluster warning path。
- 因此 source-backed boundary 是：
  - `G = 19` -> 需要 WCB 提醒
  - `G = 20` -> 不再触发 `< 20` 路径
- 一致性检查阈值来自 Python `check_clustering_consistency()`：
  `pct_clusters_with_variation < 5`，不是 `<= 5`。
- 因此 `Task 4.4` 的 exact boundary 应写成：
  - `4%` -> `is_consistent = TRUE`
  - `5%` -> `is_consistent = FALSE`

### 4. `Task 4.7` / `Task 4.8` / `Task 4.9` 的 API-scope contract

- `Task 4.7`：
  Python public `recommend_clustering_level()` 明确只保留 `ranked[1:3]`，
  即最多两个 alternatives。R 侧对应 contract 也应保持“最多 2 个”，
  这是 source-backed exact boundary。
- `Task 4.8`：
  `is_valid_cluster` / `is_valid` 的 exact 条件只有三条：
  `!is_nested_in_unit && n_clusters >= 2 && level != "lower"`。
  不应再隐含第四条“处理不得组内变化”。
- `Task 4.9`：
  必须区分两条排序语义：
  - Python internal `_generate_clustering_recommendation()`：
    双键 `(n_clusters >= 20, reliability_score)`；
  - Python public `recommend_clustering_level()`：
    单键 `reliability_score`；
  - R public `recommend_clustering()`：
    按 `.kiro/specs/story-E8-04/design.md` 保留 internal 的双键排序，
    这是明确的 `design-delta`，不是 Python bug。
- 结论：`Task 4.9` 的 source-backed regression 应优先锁定 internal helper /
  R public contract，不要误把 Python public 单键排序当成 R drift。

### 5. `Task 5.2` 的有效边界与无效边界

- 以下边界已有明确 contract，可安全写成 regression：
  - 单个候选变量；
  - 所有候选变量无效时 `recommended_var = NULL` 的 R public 设计差异；
  - `G = 1` 时 `is_valid = FALSE`；
  - `G = N` 时 `level = "same"`；
  - `CV = 0` 时 `cv_score = 1`；
  - `G_t = 0` 时 `balance_score = 0`。
- 但 `tasks.md` 草案中的 `G = 0` 目前**不是** source-backed parity truth：
  - Python `_analyze_cluster_var()` 对空 DataFrame 的最小复现结果是
    `ValueError: cannot convert float NaN to integer`；
  - 当前 R `.analyze_cluster_var()` 会返回
    `n_clusters = 0`、`min_size = 0`、`max_size = 0`、
    `reliability_score = 0.2`、`units_per_cluster = NaN` 的结构化结果。
- 论文没有给出 empty-data clustering contract，现有 spec 也未明确规定
  `G = 0` 时应报错还是返回结构化结果。
- 因此 `G = 0` 现阶段应被标记为 `undefined-boundary`：
  它足以影响后续测试设计，但证据不足以把 Python 定性为 confirmed bug。

## Bug Judgment

- `confirmed-no-new-bug`：
  本轮没有发现新的 Python-paper 冲突，不更新
  `Docs/Python包bug列表.md`。
- `suspected-but-not-ledger-worthy-yet`：
  空数据 / `G = 0` 路径目前存在 Python helper crash 与 R helper structured
  return 的行为分叉，但由于论文与 spec 都未定义该边界，暂不把它写成
  confirmed Python bug；应先由后续实现 / controller 决定目标 contract。

## Direct Implications

1. `story-worker` 可以立即为 `Task 4.1/4.3/4.4/4.7/4.8/4.9`、
   `Task 5.1/5.4` 编写 source-backed exact regressions。
2. `Task 5.2` 不应把 `G = 0` 写成 exact parity oracle；若要覆盖该路径，
   只能先做显式 contract 决策，再写测试。
3. `qa-parity` 在继续 few-cluster / hierarchical / Monte Carlo 之前，
   不需要再回头为这些 Task 4/5 边界找新的 Python oracle；关键合同已足够。

## 建议下一步

- `story-worker`：
  先补 `Task 4/5` 的 exact boundary regressions，
  并把 `G = 0` 从默认 parity 目标中剥离出来，等待单独 contract 决策。
