# E8-05 selection diagnostics source contract

## 本轮定位

- 本文件是 `story-E8-05` 的 theory-parity sidecar audit，用于给下一顺位的
  `selection_diagnostics` 实现预先冻结 source contract。
- 该文件**不**推翻当前控制面的判断：截至 `2026-03-24 15:12 CST`，
  `story-E8-04` 仍因 `AC-17` 的 package-hardening 问题占据 active slot；
  `E8-05` 只是下一顺位。

## 本轮结论

- `confirmed`: `Docs/lw2025.md` Section 4.4 只固定 selection mechanism 的理论边界，
  而不是 Python 诊断模块的具体打分启发式。论文真值是：
  1. unbalanced panel 下，demean / detrend 直接基于已观测的 pre-treatment 数据；
  2. demean 至少需要 1 个已观测 pre-treatment period；
  3. detrend 至少需要 2 个已观测 pre-treatment periods；
  4. 选择可以依赖 time-constant heterogeneity，但不能系统性依赖
     `Y_it(infinity)` 的时变 shock。
- `confirmed`: Python `selection_diagnostics.py` 有一个必须保留的语义分层：
  `missing-rate / pattern` family 使用 full-panel grid + `Y` 是否为 `NA`；
  `balance / attrition / unit_stats` family 使用 row-presence / observed-period
  语义，而不是复用 full-grid 的 `Y`-missing denominator。
- `confirmed`: Python `_is_never_treated()` 会把 `NA`、`never_treated_values`
  中的值、以及任意 `±Inf` 都视为 never-treated；R 侧不能只特判 `Inf` 而漏掉
  `-Inf`。
- `confirmed`: Python public `diagnose_selection_mechanism()` /
  `get_unit_missing_stats()` 都没有 `d` 参数。若 R 为 common-timing convenience
  保留 `d = NULL`，该接口属于 R-side extension，不构成 Python parity drift。
- `confirmed`: Python Factor 4 的 `balance_ratio` 边界是严格
  `> 0.8 -> 0`、`> 0.5 -> 10`、否则 `20`；因此 exact boundary
  `balance_ratio = 0.8` 应计 10 分，`balance_ratio = 0.5` 应计 20 分。
- `decision`: 本轮未发现新的 Python-paper 冲突，不更新
  `Docs/Python包bug列表.md`。

## 直接 source anchors

### 论文真值边界

- `Docs/lw2025.md:455-457`
  明确 unbalanced panels 下，demeaning / detrending 直接作用于已观测到的
  pre-treatment 数据；`dotY` 至少需要 1 个 pre-treatment 观测，
  `ddotY` 至少需要 2 个。
- `Docs/lw2025.md:457`
  明确 selection 允许依赖 unit-specific time-constant heterogeneity。
- `Docs/lw2025.md:457`
  同时明确 selection 不能系统性依赖 `Y_it(infinity)` 的 shocks。
- `Docs/lw2025.md:525`
  进一步说明 detrending 允许 level 和 trend 两类异质性都与 selection 相关。

### Python diagnostics contract

- `lwdid-py_v0.2.3/src/lwdid/selection_diagnostics.py:562-582`
  `_is_never_treated()`：
  `pd.isna(g)`、`g in never_treated_values`、以及 `np.isinf(g)` 都返回 `True`；
  因而 `-Inf` 与 `Inf` 都属于 never-treated。
- `lwdid-py_v0.2.3/src/lwdid/selection_diagnostics.py:587-675`
  `_compute_balance_statistics()` 使用 `data.groupby(ivar).size()` 统计
  `obs_per_unit`，并基于 `t < g` 的已存在行数计算
  `pct_usable_demean / pct_usable_detrend`。这里并未按 `Y` 的 `NA` 重建 full grid。
- `lwdid-py_v0.2.3/src/lwdid/selection_diagnostics.py:687-782`
  `_compute_attrition_analysis()` 同样基于行存在性与
  `first_obs / last_obs / n_periods` 计算 attrition，`dropout_before/after`
  由 `last_obs < g` 和 `last_obs < T_max` 决定。
- `lwdid-py_v0.2.3/src/lwdid/selection_diagnostics.py:792-1016`
  `_classify_missing_pattern()` 先构建 full-panel grid，再用
  `merged[y].isna()` 定义 `_missing`。MCAR / MAR / MNAR / UNKNOWN 的 taxonomy
  属于 Python-backed diagnostic heuristic，而不是论文 theorem。
- `lwdid-py_v0.2.3/src/lwdid/selection_diagnostics.py:1018-1167`
  `_assess_selection_risk()` 固定四因素加权启发式，并把 Factor 4 写成：
  `balance > 0.8 -> +0`、`balance > 0.5 -> +10`、否则 `+20`。另外，
  `dropout_before == 0` 时 Factor 3 自动不给分，因为条件要求
  `dropout_after > 0 and dropout_before > 0`。
- `lwdid-py_v0.2.3/src/lwdid/selection_diagnostics.py:1177-1238`
  `_compute_missing_rates()` 的 overall / by-period / by-cohort 缺失率均基于
  full-panel grid + `merged[y].isna()`。
- `lwdid-py_v0.2.3/src/lwdid/selection_diagnostics.py:1247-1354`
  `_compute_unit_stats()` 的 `n_observed / n_missing / first_observed / last_observed`
  基于 `observed_periods = unit_data[tvar].unique()`，即 row-presence 语义。
- `lwdid-py_v0.2.3/src/lwdid/selection_diagnostics.py:1362-1550`
  Python public API 只有 `gvar` / `controls` / `never_treated_values`，没有 `d`。

## 对下一顺位 story 的裁决

1. 若 `story-worker` 后续推进 `E8-05`，`missing-rate / pattern` family 必须保持
   full-grid + `Y`-`NA` 语义。
2. `balance / attrition / unit_stats` family 不应偷用 full-grid 缺失率分母；
   应按 Python 的 row-presence / observed-period 语义实现 `obs_per_unit`、
   `last_obs`、`n_pre_treatment` 与 unit-level usability。
3. `_is_never_treated()` / 边界测试必须把 `-Inf` 与 `Inf` 同时视为 never-treated。
4. `Factor 4` 的 exact boundary 不能写成
   ``>= 0.8 -> 0`` / ``< 0.5 -> 20``；source-backed 语义是
   `0.8 -> 10`、`0.5 -> 20`。
5. 若 R public API 保留 `d = NULL`，应明确写成 R convenience extension，而不是
   “Python 本来也这样”。

## 非 claims

- 不得把 MCAR / MAR / MNAR taxonomy、四因素 0-100 risk score、或
  `pct_usable_detrend < 90` 警告写成 `Docs/lw2025.md` 的 paper theorem。
- 不得把这份 sidecar audit 误写成 “`story-E8-05` 已成为 active story”；
  当前 active slot 仍由 `story-E8-04` 的 `AC-17` hardening blocker 占据。
