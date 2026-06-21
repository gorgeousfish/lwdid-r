# E8-04 AC-17 check-only staggered contract audit

- 时间：2026-03-24 21:10 CST
- 角色：theory-parity
- active_story：`story-E8-04`
- phase：`package-hardening`

## 本轮结论

- `confirmed`：fresh ordinary rerun
  `devtools::test(filter="transform-detrend|utils|validate")` 当前通过，因此
  `2026-03-24 18:24 CST` 曾写成 live blocker 的
  `test-transform-detrend.R` / `test-utils.R` / `test-validate.R`
  三条 contract 已不再是当前工作树上的 ordinary-session 红桶。
- `confirmed`：随后 fresh full
  `devtools::check(args = c("--no-manual"), error_on = "never")`
  仍返回 `1 ERROR / 5 WARNING / 3 NOTE`，但唯一 tests ERROR 已转移为
  check-time `staggered` contract：
  `test-staggered-integration.R:383`、
  `test-staggered-integration.R:531` 与
  `test-staggered-numerical.R:181`。
- `confirmed`：fresh ordinary rerun
  `devtools::test(filter="staggered-integration|staggered-numerical")`
  当前通过，说明上述三条并不是当前工作树上的 live theory/parity 失败，而是
  build/check 环境中的 hardening drift。
- `confirmed`：source audit 进一步确认
  `aggregate = "event_time"` 在本仓中属于
  “论文有 WATT / relative-time 理论对象 + R 把它提升为主 API 枚举值”的接口扩展。
  它不是 Python public `VALID_AGGREGATE` 的真值合同，也不是 Stata 当前公开语法的一部分。
  因此本轮不新增 Python bug ledger 条目。
- `confirmed`：deferred
  `aggregate + not_yet_treated -> never_treated` 的语义有论文 / Python 支撑，
  但 `lwdid_control_group_switch` 这个 warning class 名称本身是 R-local contract，
  不是 paper/Python 的同名真值。

## 直接证据

### 1. ordinary targeted rerun 已把旧 3-test 口径推成历史快照

本轮直接执行：

- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="transform-detrend|utils|validate", reporter="summary")'`

结果：

- `transform-detrend`、`utils`、`validate` 当前全部通过；
- 仅保留 `transform_detrend()` 小样本 warning，不再有 failure。

这说明 `18:24 CST` 曾写成 live blocker 的
`muffle_detrend` helper、aggregate constant 与 deferred warning-class
三条 ordinary-session contract 已经过期。

### 2. fresh full check 的 tests ERROR 已切到 check-time staggered bucket

本轮 fresh 执行：

- `Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual"), error_on = "never")'`

结果：

- `1 ERROR / 5 WARNING / 3 NOTE`

其中 tests ERROR 当前显示为：

- `test-staggered-integration.R:383`
  仍按“aggregation not implemented -> inference fields are NA”的旧 stub 期望断言；
- `test-staggered-integration.R:531`
  仍按“aggregate non-none emits stub message”的旧 stub message 合同断言；
- `test-staggered-numerical.R:181`
  仍按 `aggregate = "event_time"` 应报错的过期合同断言。

与此同时，fresh ordinary rerun

- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="staggered-integration|staggered-numerical", reporter="summary")'`

当前通过。这直接说明当前 tests ERROR 更像 build/check 环境中的 contract drift，
而不是当前工作树上的 live implementation/theory 失败。

### 3. `event_time` 的真值边界：论文对象存在，但主 API 枚举值是 R extension

论文真值：

- `Docs/lw2025.md:556-560` 直接定义了按 relative time 聚合的 `WATT(r)`；
- `Docs/lw2026.md:557-607` 直接定义了 cohort aggregate `tau_g` 与
  overall weighted effect `tau_omega`，但没有把它写成
  `aggregate = "event_time"` 这种 public API 枚举值。

R 当前实现：

- `lwdid-r/R/utils.R:52` 与 `lwdid-r/R/validate.R:65`
  都把 `"event_time"` 视为有效 aggregate；
- `lwdid-r/R/lwdid.R:1115-1118,1312-1318`
  当前已经显式接受 `aggregate = "event_time"` 并调用
  `aggregate_to_event_time()`；
- `lwdid-r/R/control_groups.R:79-82,114-118`
  明确注明：Python public aggregate 不含 `event_time`，R 的 `event_time`
  支持是 feature extension。

Python / Stata 对照：

- `lwdid-py_v0.2.3/src/lwdid/core.py:1901-1909`
  的 public `VALID_AGGREGATE` 只有 `none/cohort/overall`；
- 但 Python 仍保留内部 / 子模块路径
  `aggregate_to_event_time()`，
  见 `src/lwdid/staggered/aggregation.py:1637+` 与
  `src/lwdid/results.py:1184-1199`；
- `lwdid_stata/lwdid.sthlp:21-25,61-65`
  的公开语法没有 `aggregate()` 或 `event_time` 选项。

裁决：

- `event_time` 作为 relative-time / WATT 对象有论文来源；
- 但 `aggregate = "event_time"` 这个 public API 枚举值不属于跨语言真值，
  而是 R 的接口扩展；
- 因此 `test-staggered-numerical.R:181` 若继续把
  `event_time aggregate rejected` 当成真值，会把过期 stub 误写成 theory blocker。

### 4. deferred auto-switch warning 的语义与类名应分开裁决

测试与实现：

- `lwdid-r/tests/testthat/test-validate.R:1058-1071`
  断言：`aggregate="cohort"` + `control_group="not_yet_treated"` + 有 NT
  时，应发出 `lwdid_control_group_switch` warning，并返回 `"never_treated"`；
- `lwdid-r/R/validate.R:1873-1886`
  当前实现仍满足这条合同，而且还覆盖
  `aggregate="overall"` 与 `control_group="all_others"`。

真值来源：

- 论文与 R 注释支持的是真实语义：
  聚合效应需要 common NT control，见
  `Docs/lw2026.md:557-607` 与
  `lwdid-r/R/control_groups.R:71-85,111-118`；
- Python 也确实会做 auto-switch，见
  `lwdid-py_v0.2.3/src/lwdid/core.py:1696-1704`；
- 但 Python warning class 是 `DataWarning` / `UserWarning`，
  而不是 R 的 `lwdid_control_group_switch` 专名。

裁决：

- auto-switch 语义本身是 source-backed；
- `lwdid_control_group_switch` 这个类名是 R-local hardening contract；
- 因此若 check-time 对该 warning class 再次转红，应优先按 hardening/test harness
  处理，而不是登记为新的 theory drift。

## 裁决

1. `story-E8-04` 仍不能写成 closure-ready；fresh full check 继续红。
2. 但 theory-parity 这轮已确认：旧的 `transform-detrend/utils/validate`
   三条 blocker 口径已过期，当前 tests ERROR 应改写为
   check-time `staggered` bucket。
3. 其中 `event_time` 相关失败属于 R extension / stale stub contract 边界，
   不是新的 Python-paper 冲突；`Docs/Python包bug列表.md` 保持不变。
4. 下一轮若由 controller / correct-course / story-worker 继续推进，
   应优先修复：
   - `test-staggered-integration.R` 的旧 “aggregation not implemented / stub message” 合同；
   - `test-staggered-numerical.R` 对 `event_time` rejection 的过期断言；
   - 然后再 fresh rerun full `devtools::check(--no-manual)` 重定 AC-17 边界。
