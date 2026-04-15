# E8-04 AC-17 hardening boundary audit

## 本轮结论

- `confirmed`: `story-E8-04` 仍被 `AC-17` 阻塞，但当前 live blocker 已经是
  package-hardening / CRAN hygiene，而不是 clustering 理论、Layer 3/4 parity
  或 Python-paper drift。
- `confirmed`: 在当前工作树上，`man/aggregate_to_overall.Rd` 与
  `tests/testthat/test-sensitivity-pre-period.R` 都已通过 source-level parse smoke
  check，`nobs` 的 namespace wiring 也已在源文件中就位；因此
  `2026-03-24 15:12 CST` 以及 `2026-03-24 16:17 CST` 控制面里把这三项写成
  current blocker 的说法，至少已经部分过期。
- `confirmed`: 当前工作树 fresh
  `devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual"), error_on = "never")`
  仍未把 `AC-17` 转绿；当前可重复观察到的 hardening drift 以
  `man/compute_inference.Rd` 的 Rd 结构损坏为首，并伴随 hidden `.tmp_*`
  顶层文件、non-ASCII R 源码、undeclared imports / test dependencies，以及
  Rd / codoc / documentation hygiene 问题。
- `confirmed`: 本轮未发现新的 Python-paper 冲突，不更新
  `Docs/Python包bug列表.md`。

## 直接证据

### 当前工作树 source smoke checks

1. `tools::parse_Rd("lwdid-r/man/aggregate_to_overall.Rd")` 当前返回 `OK`，
   未再复现 `aggregate_to_overall.Rd` 的 live parse failure。
2. `parse(file = "lwdid-r/tests/testthat/test-sensitivity-pre-period.R")`
   当前返回 `OK`，未再复现 `unexpected end of input`。
3. `lwdid-r/NAMESPACE` 当前同时包含
   `S3method(nobs,lwdid_result)` 与 `importFrom(stats,nobs)`。
4. `lwdid-r/R/results.R` 当前定义了
   `nobs.lwdid_result <- function(object, ...) { object$nobs }`；
   配合 fresh `devtools::check()` 中 package load / namespace load 阶段通过，
   本轮未再复现 `object 'nobs' not found`。

### 当前工作树 fresh check 仍暴露的 hardening drift

1. `man/compute_inference.Rd` 在 build/check 阶段继续触发
   `unknown macro '\\item'`、`unexpected section header '\\value'`、
   `unexpected END_OF_INPUT` 等 Rd 解析警告。
2. `R CMD check` 继续报告 hidden `.tmp_*` 顶层文件。
3. `R CMD check` 继续报告 non-ASCII R 源文件：
   `R/clustering_diagnostics.R`、`R/plot_diagnostics.R`、
   `R/plot_event_study.R`、`R/sensitivity.R`、`R/staggered.R`。
4. `R CMD check` 继续报告未声明依赖：
   - R 代码：`MASS`、`gridExtra`、`patchwork`、`scales`
   - tests：`devtools`、`jsonlite`、`yaml`
5. `R CMD check` 继续报告 documentation hygiene 问题：
   - `to_latex_comparison` 缺 documentation entry
   - 多处 Rd cross-reference / `\usage` / codoc mismatch
   - `compute_inference.Rd` 与若干其他 Rd 仍有 markup 问题

## 裁决

1. 控制面可以继续把 `story-E8-04` 保持为 active，但不应再把
   `aggregate_to_overall.Rd`、namespace `nobs`、`test-sensitivity-pre-period.R`
   parse error 写成当前已确认的 live blocker。
2. 当前 `AC-17` 的理论边界已经足够清晰：剩余工作是 package-hardening，
   不是重新裁决论文真值、Python comparator 预期或 bug ledger。
3. 下一位 owner 若继续推进 `E8-04`，应优先处理 `compute_inference.Rd`
   与 CRAN hygiene 项，再决定是否重跑 closure-readiness 判定。
