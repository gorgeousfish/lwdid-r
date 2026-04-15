# story-E8-03 Layer 4 common-timing Monte Carlo 证据

## 结论

- `confirmed`: `e8_03_monte_carlo_common_timing_comparator.py` 以默认参数执行通过；
  Layer 4 已具备一条可重复运行的 common-timing Monte Carlo 验证链。
- `confirmed`: Python 与 R 在 50 次 replication 上的 Monte Carlo 汇总指标
  当前 exact / numeric parity 成立；`mean_att`、`mean_se`、`mean_bias`、
  `mean_abs_bias` 与 `coverage` 的汇总差异均在记录容差内。
- `confirmed`: Monte Carlo acceptance 当前通过；
  `abs(mean_bias) = 0.0152253881 <= 0.05`，
  `coverage = 0.96`，相对 95% 目标的 gap 为 `0.01 <= 0.02`。
- `confirmed`: 新增的
  `lwdid-r/tests/testthat/test-sensitivity-comprehensive-monte-carlo-parity.R`
  当前通过，说明 oracle 已被 R 侧回归测试消费，而不是停留在一次性脚本产物。

## 关键证据

1. 工件链路
   - 脚本：
     `_automation/test-artifacts/parity/e8_03_monte_carlo_common_timing_comparator.py`
   - 共享 fixture：
     `_automation/test-artifacts/parity/e8_03_monte_carlo_common_timing_fixture.csv`
   - JSON oracle：
     `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-monte-carlo-common-timing.json`
   - R 回归测试：
     `lwdid-r/tests/testthat/test-sensitivity-comprehensive-monte-carlo-parity.R`

2. Monte Carlo 设计
   - DGP：paper-backed common-timing Scenario 1C
   - replication 数：`50`
   - 每次 replication 的单位数：`160`
   - 时期数：`T = 6`
   - 首次处理期：`S = 4`
   - 控制变量：`x1`, `x2`
   - 目标指标：
     `lwdid_sensitivity(type="all")$transformation_comparison$demean_att`

3. 汇总指标
   - `mean_true_att = 6.6421532050`
   - `mean_att = 6.6573785931`
   - `mean_se = 0.4306102148`
   - `mean_bias = 0.0152253881`
   - `mean_abs_bias = 0.3542481140`
   - `coverage = 0.96`
   - `finite_att = TRUE`
   - `finite_se = TRUE`

4. Python -> R Monte Carlo parity
   - `mean_att` 汇总差异：`5.33e-15`
   - `mean_se` 汇总差异：`2.78e-16`
   - `mean_bias` 汇总差异：`5.84e-15`
   - `mean_abs_bias` 汇总差异：`6.66e-16`
   - `coverage` 汇总差异：`0`

5. warning 画像
   - comprehensive sensitivity 的 estimator comparison 分支在部分 replication
     上会触发 overlap / weight-CV warning。
   - 当前这批 warning 没有导致 replication error，也没有把目标指标
     `demean_att` / `demean_se` 推到非有限值。
   - 因此本轮将其记录为 controls / overlap 审计的后续线索，而不是 Layer 4
     common-timing Monte Carlo 的 blocker。

## 范围边界

- 本证据清除了 Layer 4 在 common-timing Scenario 1C 上“尚无可执行工件”的状态。
- 本证据没有覆盖：
  - Scenario 3C / 4C 等 misspecification Monte Carlo
  - staggered Monte Carlo
  - controls 分支 comparator 审计
  - Stata 精度 / FATAL 继承的 story-level收口

## 下一步建议

1. 复用相同 `fixture + JSON oracle + testthat` 模式，把 Layer 4 扩展到
   Scenario 3C / 4C 或 staggered DGP。
2. 优先推进 controls comparator：
   `castle` 使用 paper-valid unit-level `X_i`，
   `smoking` 先冻结 shared `X_i` fixture。
3. 若要收 Task 8.4 / 8.5 / 8.6，可直接基于这条 Monte Carlo 基线再加
   Stata 对照、FATAL 继承检查与更细的数值合理性断言。
