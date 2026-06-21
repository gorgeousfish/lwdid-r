# story-E8-03 smoking Stata 精度 comparator 证据

## 结论

- `confirmed`: qa-parity 本轮已为 `smoking` common-timing comprehensive transformation 分支补齐 Stata-backed executable oracle。
- `confirmed`: 当前 `lwdid_sensitivity(type="all")` 的 `demean/detrend` ATT 与 SE 均落在 ATT `< 1e-06`、SE `< 0.0001` 容差内。
- `bounded`: 本证据覆盖 Task `8.4` 在 `smoking` common-timing / no-controls transformation path` 的底层 `lwdid()` 精度；`8.5 / 8.6` 与更广的 hardening matrix 仍待补齐。

## 关键证据

1. comparator
   - 脚本：`_automation/test-artifacts/parity/e8_03_smoking_stata_precision_comparator.py`
   - JSON：`_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-precision.json`
   - 回归测试：`lwdid-r/tests/testthat/test-sensitivity-comprehensive-stata-precision-parity.R`

2. Stata 参考值
   - demean ATT：`-0.42217461502012654`
   - demean SE：`0.12079952386677335`
   - detrend ATT：`-0.22698869955616763`
   - detrend SE：`0.094068893843888307`

3. R vs Stata 差值
   - demean ATT abs diff：`1.17e-15`
   - demean SE abs diff：`1.8e-16`
   - detrend ATT abs diff：`2.22e-16`
   - detrend SE abs diff：`5.55e-17`
   - difference abs diff：`9.44e-16`
   - rel_diff abs diff：`9.44e-16`

## 范围边界

- 本证据不覆盖 `pre_period` / `no_anticipation` 分支的 Stata 手工重建。
- 本证据不覆盖 staggered / FATAL 继承 / 全字段数值合理性。

## 下一步建议

1. 继续把 `8.5` 做成可执行 FATAL inheritance regression，而不是停留在口头 blocker。
2. 为 `8.6` 增加全字段 finite / range audit，优先复用现有 real-data 与 Monte Carlo oracle。
