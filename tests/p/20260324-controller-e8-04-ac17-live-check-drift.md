# E8-04 AC-17 live check drift audit

## 时间

- `2026-03-24 17:28:52 CST`

## 本轮验证

- fresh rerun
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics|clustering-diagnostics-parity", reporter="summary")'`
  当前通过。
- fresh source smoke checks：
  - `tools::parse_Rd("lwdid-r/man/aggregate_to_overall.Rd")` → `OK`
  - `parse(file = "lwdid-r/tests/testthat/test-sensitivity-pre-period.R")` → `OK`
  - `grep -Fx 'S3method(stats::nobs,lwdid_result)' lwdid-r/NAMESPACE` 与
    `grep -Fx 'importFrom(stats,nobs)' lwdid-r/NAMESPACE` → `OK`
- fresh rerun
  `Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual"), error_on = "never")'`
  返回 `1 error / 8 warnings / 3 notes`。

## 关键事实

1. `story-E8-04` 的 Layer 1-4 parity 仍已收口；当前 active blocker 继续是
   `AC-17` 对应的 package-hardening，而不是 clustering 业务或 Layer 3/4 parity。
2. 旧的 live blocker 已在当前工作树过期：
   - `aggregate_to_overall.Rd` parse failure 未复现；
   - `test-sensitivity-pre-period.R` parse error 未复现；
   - namespace `nobs` wiring 已在 `NAMESPACE` / `R/results.R` 中就位。
3. 当前 fresh check 的 ERROR bucket 已收敛到 `visualization-export` suite。
   结合同窗 `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer5-release-regression.md`
   的 targeted release probe，当前 live release blocker 已进一步隔离为：
   - `TC-10.6.18`: by-cohort export literal message contract 未稳定匹配
   - `TC-10.6.34`: RI summary 未打印 `valid=999/999`
4. 当前 WARNING/NOTE 主要集中在：
   - non-ASCII R 源码：
     `R/clustering_diagnostics.R`、`R/plot_diagnostics.R`、
     `R/plot_event_study.R`、`R/sensitivity.R`、`R/staggered.R`
   - undeclared imports / test dependencies：
     `MASS`、`gridExtra`、`patchwork`、`scales`、`devtools`、`jsonlite`、`yaml`
   - Rd / documentation hygiene：
     Rd markup、cross-reference、missing doc entry `to_latex_comparison`、
     codoc mismatch、Rd `\usage` 参数不一致
   - NOTES：
     unavailable suggests / vignette dependency、visible globals、future timestamp

## 裁决

1. `story-E8-04` 继续保持 active；在 fresh `devtools::check()` 转绿前，不得切到
   `story-E8-05`。
2. 控制面必须把旧的 parse / namespace / hidden-file 叙述降级为历史快照。
3. 下一位 hardening owner 应优先处理：
   - `visualization-export` suite 的 `TC-10.6.18`
   - `visualization-export` suite 的 `TC-10.6.34`
   - non-ASCII R 源码
   - undeclared imports / test dependencies
   - Rd / codoc / documentation hygiene
