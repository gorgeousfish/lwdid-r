# story-E8-03 综合敏感性入口结构漂移修复证据

## 结论

- `resolved`: `lwdid-r/R/sensitivity.R` 中 `lwdid_sensitivity()` 与
  `print/summary/plot.lwdid_sensitivity_comprehensive` 已收敛为单一
  canonical 定义。
- `resolved`: 2093 行起的 `if (FALSE)` dead-code block 已移除，不再在
  canonical artifact 中保留多套过期入口 / S3 草稿。
- `remaining-work`: `story-E8-03` 仍未完成 parity 收口；下一步应继续推进
  Layer 3 真实数据 / 内置数据集回归与 Layer 4 Monte Carlo。

## 修复动作

1. 保留 `lwdid-r/R/sensitivity.R` 前半段已验证可工作的 comprehensive 实现。
2. 删除 2093 行起仅用于历史草稿保留的 `if (FALSE)` dead-code block。
3. 确认 `lwdid_sensitivity()` 公共接口继续显式暴露：
   - `max_anticipation`
   - `detection_threshold`

## 验证证据

### 1. source-level 结构检查

执行：

```r
e <- new.env()
sys.source("lwdid-r/R/sensitivity.R", envir = e)
names(formals(e$lwdid_sensitivity))
```

确认结果包含：

- `max_anticipation`
- `detection_threshold`

并且：

```r
b <- paste(deparse(body(e$lwdid_sensitivity)), collapse = "\n")
```

对以下模式均返回 `FALSE`：

- `lwdid_sensitivity <- function`
- `print.lwdid_sensitivity_comprehensive <- function`
- `summary.lwdid_sensitivity_comprehensive <- function`

### 2. targeted regression test

执行：

```bash
Rscript -e 'devtools::test(pkg="lwdid-r", filter="sensitivity-comprehensive")'
```

结果：

- `PASS 10`
- `FAIL 0`
- `WARN 0`

## 状态建议

- `parity-qa` 可将结构性 blocker 视为已解除。
- 下一步不应继续围绕“清理重复定义”打转，而应回到：
  - comprehensive 输出的 shared-data comparator 扩展
  - Layer 3 / Layer 4 parity 收口
