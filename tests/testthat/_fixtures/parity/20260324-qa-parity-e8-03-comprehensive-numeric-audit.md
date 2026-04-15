# story-E8-03 comprehensive 数值 hardening audit

## 结论

- `confirmed`: 三条 comprehensive real-data 路径当前未出现非有限数值。
- `confirmed`: 所有 ATT-like 标量当前均落在 `|ATT| < 10 * sd(y)` 的 hardening 边界内。
- `bounded`: 本证据覆盖 `smoking`、`castle` 与 frozen-controls 三条 `type="all"` 返回路径，主要用于收敛 Task `8.6`。

## 工件

1. hardening audit
   - 脚本：`_automation/test-artifacts/parity/e8_03_comprehensive_numeric_audit.py`
   - JSON：`_automation/test-artifacts/parity/20260324-qa-parity-e8-03-comprehensive-numeric-audit.json`
   - 测试：`lwdid-r/tests/testthat/test-sensitivity-comprehensive-hardening-parity.R`

## 案例摘要

### smoking-no-controls

- 数据路径：`smoking`
- 数值标量数：`338`
- 最大 ATT/sd(y) 比值：`1.669231`
- 非有限路径数：`0`
- 越界路径数：`0`

### castle-no-controls

- 数据路径：`castle`
- 数值标量数：`142`
- 最大 ATT/sd(y) 比值：`0.194894`
- 非有限路径数：`0`
- 越界路径数：`0`

### smoking-frozen-controls

- 数据路径：`smoking-frozen-controls`
- 数值标量数：`343`
- 最大 ATT/sd(y) 比值：`1.605923`
- 非有限路径数：`0`
- 越界路径数：`0`

## 下一步建议

1. 继续推进 Task `8.5`，把 FATAL inheritance 也固化为可执行 regression。
2. 若后续补齐 `8.4` 其余 Stata 覆盖，应把同一 hardening audit 扩展到新增 oracle 所覆盖的分支。
