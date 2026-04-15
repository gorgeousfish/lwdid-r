# story-E8-04 clustering consistency NA/scope contract

## 时间

- `2026-03-24 07:24 CST`

## Target

- 为 `story-E8-04` 的 `Task E8-04.2 / 3` 冻结 `check_clustering_consistency()`
  与 `.detect_treatment_variation_level()` 的 source-backed 真值合同，避免后续
  R 实现把 Python `nunique(dropna=True)` 语义写丢，或把 consistency check
  误写成“跨全部候选层级”的全局搜索。

## Highest-Truth Anchors

- `Docs/lw2026.md:641-643`：Lee-Wooldridge 2026 只给出“在处理变异层级聚类”
  的 cross-sectional regression 原理，并未定义 `NA` 编码或一致性函数的
  API 细节。
- `lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:697-710`：
  `_detect_treatment_variation_level()` 按 `data[v].nunique()` 升序遍历候选变量，
  返回第一个组内 treatment unique count 恒等于 `1` 的层级。
- `lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:632-640`：
  helper-level `treatment_per_cluster` 使用 `groupby(...).nunique()`，pandas
  默认 `dropna=True`。
- `lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py:1272-1280`：
  `check_clustering_consistency()` 对当前 `cluster_var` 计算
  `groupby(cluster_var)[treatment_var].nunique()`，并把
  `_detect_treatment_variation_level()` 的候选集合收窄为 `[cluster_var]`。

## Minimal Reproduction

- 复用现有 fixture：
  `_automation/test-artifacts/parity/e8_04_clustering_never_treated_fixture.csv`
- 该 fixture 的三个 cluster 的 `gvar` 取值为：
  - `A = {3}`
  - `B = {NA, 0}`
  - `C = {Inf, 4}`
- Python 实测：
  `df.groupby("cluster")["gvar"].nunique().to_dict()`
  返回 `{A: 1, B: 1, C: 2}`，因此
  `pct_with_variation = 1 / 3 * 100 = 33.333...%`。
- R 实测：
  `uniqueN(gvar, na.rm = TRUE)` 同样得到 `A=1, B=1, C=2`；
  若错误使用 `uniqueN(gvar)` 或 `length(unique(gvar))`，则会把 `B` 误算成 `2`，
  从而把不一致比例错误放大到 `2 / 3 * 100 = 66.666...%`。

## Source-Backed Contract

### 1. Unique-count semantics must stay `dropna=True`

- `.detect_treatment_variation_level()` 的候选排序和组内 treatment unique count，
  以及 `check_clustering_consistency()` 的 `cluster_treatment` 计数，都应对标
  pandas `.nunique()` 默认 `dropna=True`。
- R 侧对应实现应使用 `data.table::uniqueN(..., na.rm = TRUE)`，而不是
  `uniqueN(...)` 或 `length(unique(...))`。
- 因此在 gvar 模式下：
  - `{NA, 0}` 是 never-treated single-valued cluster，不算 variation；
  - `{Inf, 4}` 仍算 variation，因为去掉 `NA` 后仍有两个不同的非缺失值。

### 2. Consistency check is local, not global

- Python `check_clustering_consistency()` 调用
  `_detect_treatment_variation_level(data, ivar, [cluster_var], gvar, d)`，
  只对当前被检查的 `cluster_var` 做局部判断。
- 因此该函数返回的 `treatment_variation_level` 只有两种来源：
  - 若当前 `cluster_var` 组内 treatment unique count 恒为 `1`，返回 `cluster_var`
  - 否则返回 `ivar`
- 它不会扫描其他候选层级，也不会自动返回一个外部更高层级的变量名。

### 3. Paper principle and API contract must stay separated

- 论文只支持“当政策在更高层级变化时，应在该层级聚类”的识别原则。
- `check_clustering_consistency()` 的 `5%` 容差、局部 `treatment_variation_level`
  字段、以及 `dropna=True` 的计数方式，都是 Python / implementation contract，
  不是 lw2026 的显式定理。
- 后续 R 侧若要把 consistency check 改成全局搜索，需要作为显式设计变更记录，
  不能默默以“更贴近论文”为由改写 Python-backed API 语义。

## Implementation Implications

1. `Task E8-04.2` 的 `.detect_treatment_variation_level()` 应统一使用
   `data.table::uniqueN(..., na.rm = TRUE)`，包括候选排序和组内 treatment count。
2. `Task E8-04.3` 的 `check_clustering_consistency()` 必须：
   - 对 `cluster_treatment` 使用 `uniqueN(..., na.rm = TRUE)`；
   - 继续把 `.detect_treatment_variation_level()` 的候选集限制为 `c(cluster_var)`；
   - 在后续测试中用
     `_automation/test-artifacts/parity/e8_04_clustering_never_treated_fixture.csv`
     固化 `33.333...%` 而非 `66.666...%` 的 inconsistency share。
3. 该合同是 `spec-hardening`，不要求更新 `Docs/Python包bug列表.md`。

## Bug Judgment

- `confirmed-no-new-bug`: 本轮未发现新的 Python-paper 冲突。
- `confirmed`: `check_clustering_consistency()` 的 local-scope + `dropna=True`
  语义属于 Python API contract，应在 R 侧显式复现。
