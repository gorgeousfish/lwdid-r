# E8-04 AC-17 release contract scope audit

- 时间：2026-03-24 17:31 CST
- 角色：theory-parity
- active_story：story-E8-04
- phase：package-hardening

## 本轮结论

- `confirmed`：`story-E8-04` 的 Layer 1-4 parity 仍已收口；当前 fresh
  `devtools::check(args = c("--no-manual"))` 的 live ERROR 已经不再来自
  `compute_inference.Rd`、`aggregate_to_overall.Rd`、`test-sensitivity-pre-period.R`
  或 namespace `nobs`。
- `confirmed`：当前残余 ERROR 属于 release-layer contract drift，而不是论文公式、
  Python 对标或数值 parity 漂移。具体收敛为
  `visualization-export` suite 的 `TC-10.6.18` 与 `TC-10.6.34`。
- `confirmed`：`Docs/lw2025.md` / `Docs/lw2026.md` 都不定义 CSV error message 的
  literal 文案，也不规定 summary 输出必须写成 `valid=N/N` 这种 token；因此这两条
  ERROR 都应归入 package-hardening / release contract，而不是 theory blocker。
- `confirmed`：本轮未发现新的 Python-paper 冲突，不更新
  `Docs/Python包bug列表.md`。

## 直接证据

### 1. fresh check 已把旧 blocker 推成历史快照

本轮 fresh
`Rscript -e 'devtools::check(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", args = c("--no-manual"), error_on = "never")'`
当前返回 `1 ERROR / 5 WARNING / 3 NOTE`。日志中已不再出现：

- `compute_inference.Rd` parse drift
- `aggregate_to_overall.Rd` parse drift
- `test-sensitivity-pre-period.R: unexpected end of input`
- namespace `nobs` load failure

这说明 `2026-03-24 16:40 CST` 之前把这些项写成 live blocker 的结论已经再次过期。

### 2. `TC-10.6.18` 是 ASCII 会话下的 release 文案匹配问题

`lwdid-r/R/export.R` 当前在 `to_csv.lwdid_result(..., what = "by_cohort")` 的
non-staggered 分支显式抛出中文错误：

- `无队列特定结果（仅Staggered模式可用）`

在本轮 `R CMD check` 的 `session charset: ASCII` 下，targeted reproduction 得到的是：

- `<U+65E0><U+961F><U+5217><U+7279><U+5B9A><U+7ED3><U+679C><U+FF08><U+4EC5>Staggered<U+6A21><U+5F0F><U+53EF><U+7528><U+FF09>`

因此 `TC-10.6.18` 当前失败是 literal regex 无法稳定匹配 escaped message，
属于 release/encoding contract，而不是 theory 或 estimator drift。

### 3. `TC-10.6.34` 是 RI summary 格式合同未满足

用 `tests/testthat/helper-visualization-export.R` 的 common-timing mock 加上：

- `ri_pvalue = 0.004`
- `ri_distribution = seq(-3, 3, length.out = 999)`
- `ri_seed = 42`
- `ri_method = "permutation"`
- `ri_valid = 999`
- `rireps = 999`

调用 `summary()` 再 `print()`，当前输出为：

- `method=permutation | reps=999 | seed=42 | valid=999`

而不是测试要求的：

- `valid=999/999`

对应源码路径是：

- `lwdid-r/R/results.R`：`summary.lwdid_result()` 把 `object$ri_valid` 映射为
  `out$ri_n_valid`
- `lwdid-r/R/results.R`：`print.summary.lwdid_result()` 仅打印
  `sprintf("valid=%d", x$ri_n_valid)`，没有把 `rireps` 拼入同一 token

因此 `TC-10.6.34` 当前属于 release print contract drift，而不是 RI 理论错误。

### 4. 论文真值层未被当前 ERROR 推翻

`Docs/lw2025.md` 与 `Docs/lw2026.md` 当前仍只对以下内容提供真值：

- exact t inference 的理论边界
- RI p-value 的 `c/N` 定义
- transformed ATT 的识别与推断结构

它们都没有规定：

- `to_csv()` 的中文错误文案必须如何编码
- summary RI detail line 必须采用 `valid=N/N` 还是 `valid=N`

因此这两条 live ERROR 只能写成 release contract / package-hardening 问题。

## 裁决

1. controller / hardening owner 不应再把 `compute_inference.Rd`、
   `aggregate_to_overall.Rd`、`test-sensitivity-pre-period.R` 或 namespace `nobs`
   写成 current blocker。
2. 当前 `AC-17` 的 theory 结论已经稳定：剩余 ERROR 是 release-layer string /
   summary-format contract，而不是 paper/Python drift。
3. 下一轮若继续收口 `story-E8-04`，应优先修：
   - `TC-10.6.18` 的 ASCII-safe error-message contract
   - `TC-10.6.34` 的 RI summary `valid/reps` token contract
   - 其后再处理 documentation / codoc / imports / globals NOTE 桶
