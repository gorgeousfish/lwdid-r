# story-E8-03 综合敏感性入口结构漂移证据

## 结论

- `confirmed`: `lwdid-r/R/sensitivity.R` 中的综合敏感性入口与 S3 方法存在重复且互相嵌套的定义，已经形成明确的实现结构漂移。
- `suspected`: 当前活跃的 `lwdid_sensitivity()` 顶层 formals 丢失了较早实现块中的 `max_anticipation` 与 `detection_threshold` 显式参数；由于仍保留 `...`，是否构成用户可见行为回归需要由实现角色补做运行测试。
- `not-a-python-bug`: 该问题发生在 R 端，不构成 Python 相对论文的 bug 台账条目。

## 真值锚点

1. 论文与设计目标
   - `Docs/lw2025.md` Section 3 与 `Docs/lw2026.md` Section 3/8 只支持一个统一的综合敏感性入口逻辑，不支持多套互相覆盖的实现。
   - `.kiro/specs/story-E8-03/design.md` 明确将 `lwdid_sensitivity()` 设计为一个统一 wrapper，通过 `type` 分派，并通过 `...` 透传如 `max_anticipation`、`detection_threshold` 等参数。
2. Python 参考实现
   - `lwdid-py_v0.2.3/src/lwdid/sensitivity.py` 仅有一个 `sensitivity_analysis()` 主入口。
   - 同文件中的 no-anticipation 分析显式暴露 `max_anticipation` 与 `detection_threshold`。

## R 端确认到的结构漂移

### 1. 顶层重复定义

`lwdid-r/R/sensitivity.R` 中至少存在以下顶层重复定义：

- `lwdid_sensitivity <- function`：行 1365、1990、2330、2792
- `print.lwdid_sensitivity_comprehensive <- function`：行 1772、2415、3136、3284
- `summary.lwdid_sensitivity_comprehensive <- function`：行 1837、2464、3335、3477
- `plot.lwdid_sensitivity_comprehensive <- function`：行 3417

这已经偏离“单一入口 + 单一 S3 trio”的 story 设计。

### 2. 活跃定义内部再次嵌套定义

通过在独立环境中执行：

```r
e <- new.env()
sys.source("lwdid-r/R/sensitivity.R", envir = e)
deparse(body(get("lwdid_sensitivity", envir = e)))
```

可确认当前活跃的 `lwdid_sensitivity()` 函数体内部还包含：

- 局部 `.run_transformation_comparison`
- 局部 `.run_estimator_comparison`
- 局部 `print.lwdid_sensitivity_comprehensive`
- 局部 `summary.lwdid_sensitivity_comprehensive`
- 另一个嵌套的 `lwdid_sensitivity <- function(...)`

这说明文件不是“后面的定义覆盖前面的定义”这么简单，而是已经把多轮草稿拼进了当前活跃函数体。

### 3. 综合流程被实现成多套并存逻辑

活跃函数体可分出至少两套综合逻辑：

1. 外层综合流程先运行四步分析：
   - `pre_period_result` 初始化与运行
   - `no_anticipation_result` 初始化与运行
   - `transformation_comparison` 初始化与运行
   - `estimator_comparison` 初始化与运行
2. 随后函数体内部又定义一个新的 `lwdid_sensitivity <- function(...)`
3. 再后面又出现一套 transformation / estimator / overall assessment 逻辑

这意味着当前文件中并非只有“重复文本”，而是多套流程被同时保留。

### 4. 参数接口漂移为 `suspected`

活跃 `lwdid_sensitivity()` 的 formals 当前为：

```r
data, y, ivar, tvar, gvar, d, post, rolling, estimator,
controls, vce, cluster_var, type, alpha, verbose, ...
```

较早的实现块曾显式包含：

- `max_anticipation`
- `detection_threshold`

由于 story 设计允许通过 `...` 透传这些参数，这里暂记为 `suspected` 接口漂移，而不是立即定性为 confirmed 行为 bug。

## 对 parity 收口的影响

- 该问题是 R 端 `story-E8-03` 的结构性 blocker。
- 在 `story-worker` 清理为单一入口与单一 S3 方法集合之前，`qa-parity` 不应把 E8-03 视为可做最终 parity 收口。
- 当前不应更新 `Docs/Python包bug列表.md`，因为没有发现 Python 相对论文的 confirmed 冲突。

## 建议下一步

1. `story-worker`：
   - 将 `lwdid-r/R/sensitivity.R` 收敛为单一 `lwdid_sensitivity()` 顶层定义；
   - 保留单一的 `print/summary/plot.lwdid_sensitivity_comprehensive`；
   - 清理嵌套的局部 S3 定义与重复综合逻辑。
2. `qa-parity`：
   - 在结构清理后重跑 `story-E8-03` 相关测试与 comparator；
   - 重点验证 comprehensive 输出结构、recommendations 与 transformation / estimator 比较阈值逻辑。
