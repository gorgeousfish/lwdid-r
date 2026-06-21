# E8-04 closure-readiness audit

## 本轮验证

- fresh rerun
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics|clustering-diagnostics-parity", reporter="summary")'`
  当前通过。
- fresh rerun
  `Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual"), error_on = "never")'`
  当前返回 `1 error / 11 warnings / 7 notes`。

## 关键事实

- `story-E8-04` 的 Layer 3 / Layer 4 parity 当前已全部转绿，wild-bootstrap 不再是
  blocker。
- `story-E8-04` 不能被写成 closure-ready，因为 `AC-17` 依赖的 fresh
  `devtools::check()` 仍未通过。
- 当前 check 失败的主要 blocker 包括：
  - `man/aggregate_to_overall.Rd` 与 `man/compute_inference.Rd` 的 Rd 结构损坏；
  - namespace load 时报 `object 'nobs' not found`；
  - `tests/testthat/test-sensitivity-pre-period.R` 存在 parse error
    （unexpected end of input）；
  - 另有 hidden `.tmp_*` 顶层文件、non-ASCII R 源文件等 hardening 噪音。

## 结论

- 当前真实 drift 不是 clustering 业务实现，而是 control-plane 仍把
  `story-E8-04` 写成可切槽或 closure-ready。
- 在 fresh `devtools::check()` 转绿前，`story-E8-04` 仍应保持 active；`story-E8-05`
  继续保持下一顺位。
