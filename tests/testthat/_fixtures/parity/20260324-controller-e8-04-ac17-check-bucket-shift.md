# E8-04 AC-17 full-check blocker bucket shift audit

- 时间：2026-03-24 18:29 CST
- 角色：controller
- active_story：story-E8-04
- phase：package-hardening

## 本轮结论

- `confirmed`：fresh `Rscript /Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e8_04_layer5_release_regression.R` 当前通过；targeted Layer 5 release probe 已不能再单独支撑“live blocker 只剩 direct release contract”这一旧叙述。
- `confirmed`：fresh `Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual"), error_on = "never")'` 当前仍返回 `1 ERROR / 5 WARNING / 3 NOTE`；`AC-17` 仍未收口，因此 `story-E8-04` 必须继续占据 active slot，`story-E8-05` 不能抢占。
- `confirmed`：当前控制面必须把 blocker 描述从过期的 “`visualization-export` suite two-contract only” 收敛为 “full `devtools::check()` 的 broader tests/doc hardening bucket”。
- `confirmed`：`Docs/lw2025.md` / `Docs/lw2026.md` 仍只提供识别、推断与 higher-level clustering 的论文真值，不规定 CSV error literal 或 summary `valid=N/N` token；因此这些失败继续属于 package-hardening，而不是 theory drift。

## 直接证据

### 1. targeted release regression 当前通过

本轮直接执行：

- `Rscript /Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e8_04_layer5_release_regression.R`

输出为：

- `Layer 5 release regression passed.`

这说明普通 R session 下的 direct probe 已不再复现旧的 `TC-10.6.18 / TC-10.6.34` 口径。

### 2. fresh full check 仍然转红

本轮 fresh 执行：

- `Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual"), error_on = "never")'`

结果为：

- `1 ERROR / 5 WARNING / 3 NOTE`

并且 `checking tests ...` 继续报 `ERROR`。因此 `story-E8-04` 仍不能写成 closure-ready。

### 3. retained check-time test log 说明红桶明显更宽

保留在临时目录的最近 check-time `tests/testthat.Rout` 仍显示，`R CMD check` 环境下至少包含以下失败示例：

- `test-validate.R:1064`：`lwdid_data` warning class contract
- `test-vce-integration.R:109`：`HC3 full vcov matches sandwich`，当前是 dimnames contract drift
- `test-vce-return-values.R:109`：当前返回字段数 `21` vs 旧测试期望 `19`
- `test-visualization-export.R:402`：check-time `TC-10.6.18`
- `test-visualization-export.R:568`：check-time `TC-10.6.34`

与此同时，本轮 local targeted rerun

- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="vce-integration", reporter="summary")'`

当前通过。该对照说明 full-check 红桶不是一个可以被单个 local suite 代表的简单 blocker，而是 build/check 环境下的 broader contract drift。

### 4. 文档与 NOTE 桶仍在

fresh full check 仍保留：

- `R code for possible problems` note
- `Rd files`
- `Rd cross-references`
- `missing documentation entries`
- `code/documentation mismatches`
- `Rd \usage sections`

因此即便 targeted release probe 通过，`AC-17` 也仍未满足。

## 裁决

1. `phase` 保持 `package-hardening`。
2. `active_story` 保持 `story-E8-04`。
3. `current_goal` 与 parity / finding / spec 注记都不得继续把 live blocker 简化成“只剩 `TC-10.6.18 / TC-10.6.34`”。
4. 下一顺位应继续让 hardening owner 基于 fresh full-check 环境处理 tests/doc buckets；controller 本轮只收敛控制面，不替 worker 抢实现。
