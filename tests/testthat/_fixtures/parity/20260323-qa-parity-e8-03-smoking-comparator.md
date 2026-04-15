# story-E8-03 smoking 内置数据 comparator 证据

## 结论

- `confirmed`: `smoking.csv` 的 shared-data Python -> R comparator 已通过；
  `story-E8-03` 在 common-timing / no-controls 的真实内置数据路径上，
  exact / numeric parity 当前成立。
- `confirmed`: 本轮新增的
  `lwdid-r/tests/testthat/test-sensitivity-comprehensive-realdata-parity.R`
  已通过，说明 comparator oracle 已被 R 侧回归测试消费，而不是孤立留在
  `_automation/test-artifacts/parity/`。
- `resolved`: `parity-qa-status.yaml` 中关于 `lwdid-r/R/sensitivity.R`
  “仍存在多套顶层综合入口 / S3 定义”的 blocker 已过时；当前源码只剩一套
  canonical 顶层定义。

## 关键证据

1. shared-data comparator
   - 脚本：
     `_automation/test-artifacts/parity/e8_03_smoking_realdata_comparator.py`
   - 报告：
     `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-comparator.json`
   - 数据源：
     `lwdid-py_v0.2.3/data/smoking.csv`
   - 样本：1209 行、39 个州、1 个 treated unit、1970-2000

2. Python -> R 字段级一致性
   - `pre_period.is_robust = FALSE`
   - `pre_period.robustness_level = "sensitive"`
   - `anticipation.recommended_exclusion = 0`
   - `pre_period.sensitivity_ratio` 绝对差
     `1.22e-15`，容差 `1e-6`
   - `demean_att / detrend_att / difference` 全部在 `1e-6` 内
   - `demean_se / detrend_se` 全部在 `1e-4` 内
   - `estimator` 在 Python 与 R 端都为 `NULL`，符合 no-controls 路径预期

3. warning 画像
   - Python 与 R 都会在该数据上重复触发 small-sample warning，
     根因是 `N_treated = 1`
   - R 回归测试已改为内部捕获 warning，并对照 oracle 的 warning 文本，
     因此 targeted test 输出保持 `WARN 0`

4. 论文 / Stata 旁证
   - 当前 R 输出相对 Table 3 参考值的差值：
     - demean ATT：`1.746e-4`
     - demean SE：`2.005e-4`
     - detrend ATT：`8.27e-9`
     - detrend SE：`5.09e-9`
   - 这说明 smoking real-data comparator 同时为 Layer 0 论文再现提供了
     一条可复核旁证，但尚未覆盖 README/example 全面再现。

## 范围边界

- 本证据推进了 Layer 3 的 common-timing 内置数据验证。
- 本证据没有覆盖：
  - `castle` / staggered 真实数据 comparator
  - Layer 4 Monte Carlo
  - comprehensive sensitivity 在 time-varying controls 路径下的额外分叉

## 下一步建议

1. 复用同一 comparator 结构，为 `castle` 建立 staggered real-data oracle。
2. 在 Layer 4 增加 Monte Carlo acceptance 规则与可执行脚本。
3. 如需推进 comprehensive sensitivity 的 controls 路径，再单独审计
   smoking 上 time-varying controls 的 warning / fallback 行为。
