# story-E8-04 Layer 5 release regression

## 时间
- `2026-03-27 09:11:09 CST`

## 范围
- 重跑 `e8_04_package_hardening_regression.R`，确认旧的 `Rd/nobs/parse` hardening probe 是否仍失败。
- 运行 `devtools::test(filter = "visualization-export")` 作为当前 release blocker suite。
- 直接复现 `to_csv(..., what = "by_cohort")` 与 RI summary print contract。

## 关键结果
- legacy hardening probe: `passed`。
- visualization-export suite: `passed`。
- suite failure labels: `none`。
- by_cohort export literal match: `true`；actual message: `No cohort-specific results (Staggered mode only).`。
- RI summary header present: `true`；valid fraction present: `true`；detail line: `  method=permutation | reps=999 | seed=42 | valid=999/999`。

## 产出文件
- JSON: `/Users/cxy/Desktop/lwdid_r/lwdid-r/tests/testthat/_fixtures/parity/20260324-qa-parity-e8-04-layer5-release-regression.json`
- 说明: `/Users/cxy/Desktop/lwdid_r/lwdid-r/tests/testthat/_fixtures/parity/20260324-qa-parity-e8-04-layer5-release-regression.md`

## 裁决
- 旧的 `aggregate_to_overall.Rd` / `compute_inference.Rd` / `test-sensitivity-pre-period.R` / namespace `nobs` blocker 已退出 live blocker 集。
- 当前 targeted Layer 5 release probe 已转绿：`visualization-export` suite、`by_cohort` export literal contract 与 RI summary `valid=999/999` contract 当前均通过。
- 后续若 `story-E8-04` 仍未 closure-ready，剩余 blocker 应由 fresh `devtools::check(--no-manual)` 的 test / Rd / documentation hygiene 结果单独维护，不得继续把 `TC-10.6.18` / `TC-10.6.34` 写成 live blocker。
- 该证据只覆盖 targeted release blocker probe；完整 `devtools::check(--no-manual)` 仍需单独作为 closure-readiness 证据维护。
