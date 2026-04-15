# story-E8-03 controls comparator 小样本 simple-controls 审计

## 结论

- `confirmed`: `castle` 的 controls 路径不是因为 covariates 时变而失效。
  `income`、`unemployrt`、`poverty` 在 `sid` 内均为 time-invariant，符合论文中的
  unit-level `X_i` 设定。
- `confirmed`: 论文 `lw2025/lw2026` 方程 `(2.18)` 明确允许在
  `N > K + 2` 时估计 simple-controls OLS；R 端 `estimate_ra_common()`
  也按这一 paper-backed tier-2 fallback 执行。
- `confirmed`: `lwdid-py_v0.2.3` 在 controls 样本不满足
  `N_1 > K+1` / `N_0 > K+1` 时，会直接忽略 controls，而不是退到
  paper-valid 的 simple-controls 回归。`castle` 上这会把 controls 路径退回
  no-controls 数值，因此不应被当作 parity 真值。
- `boundary`: `smoking` 常见 raw covariates
  `lnincome` / `beer` / `age15to24` / `lretprice` 确实是 time-varying `X_it`。
  该数据集如要做 controls comparator，必须先冻结为 shared unit-level `X_i`
  fixture，不能直接拿 raw panel 列做 source-backed comparator。

## 真值锚点

### 1. 论文中的 covariate 对象始终是 unit-level `X_i`

- `Docs/lw2025.md:154-158` 与 `Docs/lw2026.md:160-164` 的方程 `(2.18)` 都写成
  `\Delta \bar{Y}_i = \alpha + \tau D_i + \mathbf{X}_i \beta + U_i`，
  并明确只要求 `N > K + 2` 即可做 OLS。
- `Docs/lw2025.md:477-479` 把 RA / IPWRA 的条件均值对象写成
  `\dot{m}_{rg}(\mathbf{X})`，说明这里的 controls 仍是 unit-level `X_i`，
  不是 period-specific `X_{it}`。

结论：只要 controls 是合法的 `X_i`，full interaction 不可行时仍不能自动把
 controls 整段删掉；`(2.18)` 的 simple-controls 回归仍然是 paper-backed 备选。

### 2. 两端源码对 `X_i` 合法性的共同约束

- Python `validation.py` 明确写道 controls 必须是 time-invariant，
  “requires time-constant controls `X_i`, not time-varying `X_it`”。
- R `validate.R` 也对 controls 做同样的 time-invariant 检查。
- R `staggered.R` 与 Python `staggered/estimation.py` 的接口注释都把
  staggered controls 描述为 time-invariant control columns。

因此，`smoking` raw controls 的问题是“输入不合法”；`castle` 这组 controls
 则是合法 `X_i`，应进入 paper-backed 小样本回归比较。

## 源码证据

### R 当前行为

- `lwdid-r/R/estimate.R`
  - Tier 1：当 `N_1 > K+1` 且 `N_0 > K+1` 时使用 full interaction；
  - Tier 2：当总样本满足 `N > K + 2` 时，降级为 simple controls，
    注释直接标为 `lw2026 eq. 2.18`；
  - Tier 3：只有在 `N <= K + 2` 时才丢弃所有 controls。

### Python 当前行为

- `lwdid-py_v0.2.3/src/lwdid/estimation.py`
  的 `_prepare_controls_for_ra()` 只在
  `N_treated > K+1` 且 `N_control > K+1` 时纳入 controls。
- 否则它直接发出
  `Controls will be ignored in the regression.`
  的 warning，并返回 `include = False`，没有 paper-backed 的 simple-controls
  fallback。

## 最小复现

### A. `castle` controls 路径：合法 `X_i` 上的 Python/R 漂移

数据与 controls：

- 数据：`lwdid-py_v0.2.3/data/castle.csv`
- `gvar = effyear`
- controls：`income`, `unemployrt`, `poverty`
- 校验：三者在 `sid` 内的最大 within-unit std 都为 `0.0`

fresh 复跑 `lwdid_sensitivity(type = "all", verbose = FALSE)` 后：

- Python
  - warnings 中出现
    `Controls not included: N_treated=1, N_control=29, K=3. Need N_treated > K+1 and N_control > K+1.`
  - `estimator_comparison = {ra = 0.091745..., ipw = 0.091745..., ipwra = 0.091745..., range = 0.0}`
  - `overall_assessment = "Multiple concerns: pre-period sensitivity, anticipation effects, transformation sensitivity..."`
- R
  - warnings 中多次出现
    `degraded to simple controls (lw2026 eq. 2.18)`
  - `estimator_comparison = {ra = 0.0645628..., ipw = 0.0645628..., ipwra = 0.0645628..., range = 0, rel_range = 0, baseline_att = 0.0645628...}`
  - `overall_assessment = "多项问题：前期敏感性、预期效应，需谨慎解读"`

裁决：

- Python 的 `0.091745...` 正好等于既有 no-controls comparator 的 demean ATT，
  说明该路径实质上已退回 no-controls。
- R 的 `0.0645628...` 来自 simple-controls fallback，且该 fallback 与方程
  `(2.18)` 一致。
- 因此这里不是“R 偏离 Python”，而是 Python 偏离论文允许的 controls 回归层级。

### B. `smoking` raw controls：时变 `X_it` 仅能作为边界告警，不是 parity 真值

数据与 raw controls：

- 数据：`lwdid-py_v0.2.3/data/smoking.csv`
- controls：`lnincome`, `beer`, `age15to24`, `lretprice`

在 `state` 内这些列都明显时变；fresh 复跑时：

- Python 会反复报出
  `Control variables must be time-invariant ... X_i, not time-varying X_it`
- R 侧综合敏感性各 specification 也会因同一原因失败

裁决：

- `smoking` raw controls 不是合格的 source-backed comparator 输入。
- 若后续要补 common-timing controls comparator，必须先构造 shared frozen
  controls 列，再同时喂给 Python 与 R。

## 裁决

- `Python bug status`: `confirmed`
- `bug 类型`: 小样本 controls 路径错误跳过 paper-valid simple-controls regression
- `影响范围`:
  - `story-E8-03` 的 `castle` controls / estimator-comparison parity；
  - 所有使用合法 unit-level controls、但 `N_1 <= K+1` 或 `N_0 <= K+1`、
    同时 `N > K+2` 的 Python RA / sensitivity 路径。
- `R 是否跟随`: 否。R 应继续保留 `eq. 2.18` 的 simple-controls fallback。

## 对后续 comparator 的要求

1. `castle` controls comparator：
   继续使用 unit-level `income` / `unemployrt` / `poverty`，但不得把 Python
   no-controls fallback 当真值；需先参照 bug ledger 标记 waived/drift。
2. `smoking` controls comparator：
   先明确定义 frozen `X_i` 生成规则，再做 shared-data parity。
3. `story-E8-03` 的 8.4 / controls audit：
   需要显式区分
   - paper-backed simple-controls path；
   - Python bug path；
   - truly invalid raw `X_it` path。
