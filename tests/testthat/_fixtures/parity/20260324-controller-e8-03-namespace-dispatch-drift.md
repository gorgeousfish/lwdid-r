# story-E8-03 namespace / S3 dispatch drift

## 时间

- `2026-03-24 04:11 CST`

## 复核命令

- `Rscript -e 'library(lwdid); cat("exported:", "lwdid_sensitivity" %in% getNamespaceExports("lwdid"), "\n")'`
- `Rscript -e 'pkgload::load_all("/Users/cxy/Desktop/lwdid_r/lwdid-r", quiet=TRUE); data("smoking", package="lwdid"); res <- lwdid_sensitivity(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", type="all", verbose=FALSE); cat("print_head:", paste(capture.output(print(res))[1:2], collapse=" | "), "\n"); cat("summary_head:", paste(capture.output(summary(res))[1:2], collapse=" | "), "\n")'`
- `rg -n "S3method\\(print|S3method\\(summary|export\\(lwdid_sensitivity" /Users/cxy/Desktop/lwdid_r/lwdid-r/NAMESPACE`

## 关键事实

1. `lwdid-r/R/sensitivity.R` 已定义 `lwdid_sensitivity()`、`print.lwdid_sensitivity_comprehensive()`、
   `summary.lwdid_sensitivity_comprehensive()`，且三者均带 `@export`。
2. `lwdid-r/NAMESPACE` 当前仍缺 `export(lwdid_sensitivity)`，也缺
   `S3method(print,lwdid_sensitivity_comprehensive)` 与
   `S3method(summary,lwdid_sensitivity_comprehensive)`。
3. 安装态命名空间检查返回 `exported: FALSE`，说明 `library(lwdid)` 下
   `lwdid_sensitivity()` 仍不可见。
4. 在 `pkgload::load_all()` 环境下，generic `print(res)` 的头两行仍为
   `$pre_period_result | Pre-treatment Period Robustness Analysis`，而
   `summary(res)` 的头两行仍为默认 `Length / Class / Mode` 结构；这表明 generic
   dispatch 仍退回默认 list 输出，而不是 comprehensive 的定制格式。

## 裁决

- `confirmed`: Task `E8-03.1` 仍存在 public API export drift，不能视为 closure-ready。
- `confirmed`: `AC-23` / `AC-24` 继续保持未勾选，Task `E8-03.2` 应回退为 partial。
- `confirmed`: `story-E8-03` 仍应保持 active；当前直接缺口不再是“只剩 Task 8.4”，而是
  Task `8.4` 的 `pre_period n_pre=1` Stata `r(2001)` 边界，加上本文件记录的
  namespace / S3 dispatch drift。
