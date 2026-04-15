# E8-04 hardening static convergence

## 本轮验证

- fresh rerun
  `Rscript /Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e8_04_package_hardening_static_regression.R`
  当前通过。
- fresh rerun
  `Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual", "--no-tests", "--no-examples"), error_on = "never", document = FALSE)'`
  当前返回 `5 warnings / 3 notes / 0 errors`。

## 已解除的 live drift

- `checking code files for non-ASCII characters` 当前为 `OK`。
- `checking dependencies in R code` 当前为 `OK`。
- `checking for unstated dependencies in 'tests'` 当前为 `OK`。

## 当前剩余 warning

- `Rd files`
- `Rd cross-references`
- `missing documentation entries`
- `code/documentation mismatches`
- `Rd \usage sections`

## 结论

- non-ASCII R 源码、undeclared imports 与 tests dependencies 当前不应再被写成
  `story-E8-04` 的 live blocker。
- `story-E8-04` 仍停留在 `package-hardening`，且仍不能标记为 closure-ready。
- 本轮未发现新的 Python-paper 冲突，`Docs/Python包bug列表.md` 保持不变。
