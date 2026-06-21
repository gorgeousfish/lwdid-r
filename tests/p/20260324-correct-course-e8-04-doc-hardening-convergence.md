# 2026-03-24 correct-course E8-04 文档 hardening 收敛审计

- 时间：`2026-03-24 21:10 CST`
- phase：`package-hardening`
- active_story：`story-E8-04`

## 本轮命令

1. `Rscript -e 'devtools::load_all("/Users/cxy/Desktop/lwdid_r/lwdid-r", quiet = TRUE); testthat::test_file("/Users/cxy/Desktop/lwdid_r/lwdid-r/tests/testthat/test-sensitivity-comprehensive.R", reporter = testthat::SummaryReporter$new())'`
2. `Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual", "--no-tests", "--no-examples"), document = FALSE, error_on = "never")'`

## 关键观察

- 先在 `test-sensitivity-comprehensive.R` 新增 regression，锁定
  `.new_lwdid_sensitivity_comprehensive(overall_assessment = NULL, recommendations = NULL)`
  仍应回落到默认中文总结；新增断言在修复前为红，修复后 fresh rerun 转绿。
- 随后对 canonical roxygen 源做最小 hardening：
  - 修正 `R/export.R`、`R/sensitivity.R`、`R/transform.R`、
    `R/aggregate.R`、`R/propensity.R` 中的 brace/link 误触发；
  - 为 `to_latex_comparison()` 补齐导出文档；
  - 为 `lwdid()`、`new_lwdid_result()` 与 `lwdid-conditions` 收敛参数文档；
  - 移除 `%||%` 的独立 Rd，避免 `grapes-or-or-grapes.Rd` 的 checkRd name warning；
  - 将 `.new_lwdid_sensitivity_comprehensive()` 的默认文案收口到函数体内部，
    消除 ASCII check 下的 codoc default-value 漂移。
- `devtools::document()` 后，fresh partial
  `devtools::check(--no-manual, --no-tests, --no-examples, document = FALSE)`
  当前返回 `0 ERROR / 0 WARNING / 3 NOTE`。

## 裁决

- `Rd files`
- `Rd cross-references`
- `missing documentation entries`
- `code/documentation mismatches`
- `Rd \usage sections`

以上五个文档 hardening warning bucket 当前已全部转绿，不应再被写成 `story-E8-04`
的 live blocker。

- 当前 remaining bucket 仅确定保留：
  - `R code for possible problems` NOTE
  - 外部仓库 / future timestamp 两条环境 NOTE
- 本轮没有 fresh 重放完整 `devtools::check(--no-manual)` 的 tests 阶段，因此
  `AC-17` 仍不能写成完成；`story-E8-04` 继续保持 active。
