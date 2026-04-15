# story-E8-03 castle 内置数据 comparator 证据

## 结论

- `confirmed`: `castle.csv` 的 shared-data Python -> R comparator 已通过；
  `story-E8-03` 在 staggered / no-controls 的真实内置数据路径上，
  exact / numeric parity 当前成立。
- `confirmed`: 本轮新增的
  `lwdid-r/tests/testthat/test-sensitivity-comprehensive-staggered-realdata-parity.R`
  已通过，说明 `castle` oracle 已被 R 侧回归测试消费，而不是停留在一次性 JSON。
- `confirmed`: `castle` 的 demean / detrend overall ATT 与论文 Table 3 旁证
  仅相差 `2.55e-4` 与 `4.50e-4`，可作为 Layer 0 论文再现的附属数值支撑。

## 关键证据

1. shared-data comparator
   - 脚本：
     `_automation/test-artifacts/parity/e8_03_castle_realdata_comparator.py`
   - 报告：
     `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-castle-comparator.json`
   - 数据源：
     `lwdid-py_v0.2.3/data/castle.csv`
   - 样本：550 行、50 个州、21 个 treated units、29 个 never-treated units、
     5 个 adoption cohorts（2005-2009）

2. Python -> R 字段级一致性
   - `pre_period.is_robust = FALSE`
   - `pre_period.robustness_level = "sensitive"`
   - `anticipation.anticipation_detected = TRUE`
   - `anticipation.recommended_exclusion = 2`
   - `pre_period.sensitivity_ratio` 绝对差
     `1.11e-15`，容差 `1e-6`
   - `demean_att / detrend_att / difference` 全部在 `1e-6` 内
   - `demean_se / detrend_se` 全部在 `1e-4` 内
   - `estimator` 在 Python 与 R 端都为 `NULL`，符合 no-controls 路径预期

3. warning 画像
   - Python 与 R 都会在逐步排除前期时暴露 unbalanced-panel 与
     `not_yet_treated -> never_treated` 自动切换诊断
   - 两端 warning 文案粒度不同：Python 额外展开 selection-mechanism 解释，
     R 侧则压缩为本地诊断消息
   - R 回归测试当前消费的是 oracle 中的 R warning 画像，确保 targeted test
     对当前包内 contract 稳定

4. 论文 / README 旁证
   - 当前 R 输出相对 Table 3 参考值的差值：
     - demean ATT：`2.546e-4`
     - detrend ATT：`4.497e-4`
   - 这说明 `castle` comparator 同时为 Layer 0 提供了一条可复核的 paper-side
     numerical sanity check

## 范围边界

- 本证据推进并基本补齐了 Layer 3 的 `castle` / staggered real-data 验证。
- 本证据没有覆盖：
  - Layer 4 Monte Carlo
  - `plot-diagnostics` 当前 gate 暴露的可视化回归
  - `castle` controls / estimator-comparison 分支的单独 parity 审计

## 下一步建议

1. 保持 `real_data_status` 的新基线，把 qa-parity 主力切到 Layer 4 Monte Carlo。
2. 若要继续扩大 `castle` 覆盖面，再单独审计 controls 路径，避免把 estimator
   支持差异混进当前 no-controls comparator 结论。
3. 与 `story-worker` / `controller` 协同，把 `plot-diagnostics` 的 source-backed
   blocker 清零。
