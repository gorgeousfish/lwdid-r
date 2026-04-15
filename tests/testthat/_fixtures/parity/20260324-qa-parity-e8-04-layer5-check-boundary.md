# story-E8-04 Layer 5 full-check boundary

## 时间
- `2026-03-24 18:24:56 CST`

## 范围
- 复跑 `_automation/test-artifacts/parity/e8_04_layer5_release_regression.R`，确认 targeted release probe 当前继续通过。
- fresh rerun `devtools::test(filter = "vce-return-values")` 与 `devtools::test(filter = "vce-integration")`，验证 VCE 相关 live test blocker 是否已解除。
- fresh rerun `devtools::check(args = c("--no-manual"), error_on = "never")`，重定 Layer 5 的 live ERROR/WARNING/NOTE 边界。

## 关键结果
- Layer 5 targeted release probe：`passed`。
- `vce-return-values`：`passed`；`T5e-01e` 已按 `estimate_ra_common()` 当前 21 字段返回值合同转绿。
- `vce-integration`：`passed`；`T5-02d` 已按 `lm(y ~ D)` reference 对齐 HC3 vcov 的名字合同。
- fresh `devtools::check(--no-manual)`：`1 ERROR / 5 WARNING / 3 NOTE`。

## 新的 live blocker
- `TC-10.6.18` / `TC-10.6.34` 已退出 live blocker 集；Layer 5 release contracts 当前转绿。
- VCE stale tests 也已退出 live blocker 集；当前不应再把 `vce-return-values` / `vce-integration` 写成 live ERROR。
- 当前 fresh check 的 tests ERROR 已收敛为三条：
  - `test-transform-detrend.R:378`：`muffle_detrend` helper 缺失。
  - `test-utils.R:37`：`LWDID_VALID_AGGREGATE_LEVELS` 仍按旧合同断言，未纳入 `event_time`。
  - `test-validate.R:1064`：`aggregate + not_yet_treated` 的 deferred path 不再发出预期的 `lwdid_data` warning。

## 仍存 hygiene 桶
- WARNING：`Rd files`、`Rd cross-references`、`missing documentation entries`、`code/documentation mismatches`、`Rd usage sections`。
- NOTE：`R code for possible problems`、`package dependencies`、`future file timestamps`。

## 修改文件
- `lwdid-r/tests/testthat/test-vce-return-values.R`
- `lwdid-r/tests/testthat/test-vce-integration.R`
- `_automation/test-artifacts/parity/e8_04_layer5_release_regression.R`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer5-release-regression.md`

## 裁决
- Layer 5 当前应视为 `in-progress`，而不是 `pending`。
- `story-E8-04` 仍不能写成 closure-ready，但当前 controller / correct-course / qa-parity 都不得继续把 release contracts 或 VCE stale tests 记为 live blocker。
- 下一轮若继续在 package-hardening 轨道上推进，应优先在 `transform-detrend` helper、aggregate 常量合同与 deferred validation warning contract 三处做 source audit + targeted rerun。
