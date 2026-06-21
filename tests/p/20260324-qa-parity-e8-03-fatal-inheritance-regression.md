# story-E8-03 FATAL inheritance regression

## 时间

- `2026-03-24 03:11:05 CST`

## 范围

- 使用 `generate_staggered_panel(seed = 42)` 对 `story-E8-03` comprehensive sensitivity wrapper 做 Layer 5 FATAL inheritance regression。
- 本轮覆盖 `FATAL-001` 的 `not_yet_treated` strict-mask 数值路径，以及 `FATAL-002` / `FATAL-004` 的 aggregation auto-switch + NT weighted-average 路径。
- `FATAL-003` 仍未形成 story-level `ri=TRUE` executable regression，本证据只把 blocker 缩窄到该项。

## 关键结果

- strict-mask baseline：wrapper = `3.4980540427941`；direct NYT = `3.4980540427941`；direct NT = `3.65014850962517`。
- strict-mask 判定：`wrapper_matches_direct_nyt = true`；`wrapper_differs_from_direct_nt = true`。
- aggregation 判定：`switch_warning_detected = true`；`wrapper_matches_direct_nt = true`。
- overall/demean：wrapper ATT/SE = `3.68161107583395` / `0.706206824323679`。
- overall/detrend：wrapper ATT/SE = `2.76135634921692` / `0.469865777033168`。

## 产出文件

- JSON oracle: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260324-qa-parity-e8-03-fatal-inheritance-regression.json`
- 当前说明: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260324-qa-parity-e8-03-fatal-inheritance-regression.md`

## 裁决

- `exact_status = passed`
- `numeric_status = passed`
- `remaining_gap = FATAL-003 story-level RI regression`
