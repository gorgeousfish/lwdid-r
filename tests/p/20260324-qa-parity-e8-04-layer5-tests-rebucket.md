# story-E8-04 Layer 5 tests rebucket

## 时间
- `2026-03-24 21:04:05 CST`

## 范围
- 修正 `test-transform-detrend.R`、`test-utils.R`、`test-validate.R` 三条 stale Layer 5 test contract。
- fresh rerun `devtools::test(filter = "transform-detrend|utils|validate")` 与 `e8_04_layer5_release_regression.R`。
- fresh rerun `devtools::check(args = c("--no-manual"), error_on = "never")`，重定 canonical full-check blocker。
- 对照 rerun `devtools::test(filter = "staggered-integration")` 与 `devtools::test(filter = "staggered-numerical")`，判断新 surfaced blocker 是否只在 check-time 复现。

## 关键结果
- 旧的三条 live tests blocker 已转绿：
  - `test-transform-detrend.R:378` 不再因 `muffle_detrend()` helper 缺失失败。
  - `test-utils.R:37` 已按当前 `LWDID_VALID_AGGREGATE_LEVELS = c("none", "cohort", "overall", "event_time")` 对齐。
  - `test-validate.R:1064` 已按真实 warning class `lwdid_control_group_switch` 对齐。
- `Rscript _automation/test-artifacts/parity/e8_04_layer5_release_regression.R` 当前继续通过。
- fresh `devtools::check(--no-manual)` 仍返回 `1 ERROR / 5 WARNING / 3 NOTE`，但 live tests bucket 已从旧的 `transform/utils/validate` 前移为新的 3 条 check-time failure：
  - `test-staggered-integration.R:383`：aggregate non-none path 不再保持 `ci_upper` 等统计量为 `NA`。
  - `test-staggered-integration.R:531`：aggregate non-none path 未再发出旧的 stub message contract。
  - `test-staggered-numerical.R:181`：`event_time` aggregate 不再触发旧的 reject/error contract。
- 同窗 fresh `devtools::test(filter = "staggered-integration")` 与 `devtools::test(filter = "staggered-numerical")` 当前都通过，说明这三条是 `R CMD check` 安装环境下才会浮现的 hardening drift，而不是普通 source-test 会话直接红。

## 仍存 hygiene 桶
- WARNING：`Rd files`、`Rd cross-references`、`missing documentation entries`、`code/documentation mismatches`、`Rd usage sections`。
- NOTE：`package dependencies`、`future file timestamps`、`R code for possible problems`。

## 裁决
- `layer_5` 继续保持 `in-progress`。
- `story-E8-04` 仍不得写成 closure-ready，也不得切换到 `story-E8-05`。
- canonical blocker 边界现应改写为“check-time staggered aggregation contracts + globals/Rd hygiene”，不得再把 `transform-detrend` / `utils` / `validate` 三条 stale contract 写回 live blocker。
