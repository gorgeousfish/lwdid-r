# E8-04 AC-17 check-time tests bucket drift（correct-course）

## 时间

- `2026-03-24 21:20 CST`

## 本轮执行

1. fresh rerun ordinary targeted tests：
   `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="staggered-integration|staggered-numerical|conditions", reporter="summary")'`
2. fresh rerun full hardening gate：
   `Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual"), error_on = "never")'`
3. 复核该次 full check 生成的 `00check.log` 与 `tests/testthat.Rout`

## 关键事实

- ordinary `staggered-integration|staggered-numerical|conditions` 当前通过；本轮已修正 3 条旧测试合同：
  - 默认 staggered 已返回 cohort-level inference，而非 `NA` stub
  - `aggregate="cohort"` 不再发 stub message
  - `aggregate="event_time"` 当前正常返回 `event_time_effects`
- 条件系统已把历史 `lwdid_invalid_param` 纳入 `lwdid_invalid_parameter` 继承链，因此旧代码路径不再与公共 condition 体系脱节
- fresh full `devtools::check(--no-manual)` 在 20 分钟超时前，`00check.log` 已显示：
  - `checking R code for possible problems ... OK`
  - `checking Rd files ... OK`
  - `checking Rd cross-references ... OK`
  - `checking for missing documentation entries ... OK`
  - `checking for code/documentation mismatches ... OK`
  - `checking Rd \usage sections ... OK`
- 同一次 full check 的 `tests/testthat.Rout` 已出现持续保存 `_problems/...`：
  - `_problems/test-clustering-diagnostics-parity-*`
  - `_problems/test-export-*`
- 因此当前 live blocker 已不能继续写成：
  - 旧的 `transform-detrend|utils|validate` 三条 stale contract
  - 仅剩 `Rd/doc` hardening bucket
  - 稳定收敛为 3 条 ordinary `staggered` 失败

## 结论

- `story-E8-04` 仍不得 closure-ready
- `AC-17` 的 canonical blocker 现应收敛为：
  `installed-package / full-check tests bucket`
- 现有证据只足以确认：
  - pre-tests 的 code/Rd/doc hardening 桶当前已转绿
  - full-check tests 阶段仍存在更广的 check-time failures
- 由于本轮 full check 被超时截断，当前还不能把 tests bucket 过早压缩成最终枚举列表

## 对 Python / theory 的裁决

- 本轮不更新 `Docs/Python包bug列表.md`
- 依据：
  - `aggregate="event_time"` 属于 R extension public API，不是 Python public `VALID_AGGREGATE` 真值
  - ordinary suite 已通过，当前问题发生在 build/check hardening 场景，而非 paper/Python 数值或语义漂移
