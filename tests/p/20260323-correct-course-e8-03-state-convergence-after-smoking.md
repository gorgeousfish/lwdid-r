# story-E8-03 状态收敛复核（smoking parity 之后）

## 结论

- `resolved`: `_automation` 控制平面继续把 `sensitivity.R` 结构漂移当作 blocker 的表述已经过时。
- `confirmed`: 当前源码层面只剩一套 canonical `lwdid_sensitivity()` 与 comprehensive S3 trio。
- `confirmed`: `smoking` real-data parity oracle 已存在，`test-sensitivity-comprehensive-realdata-parity` 在本轮 fresh run 中 `FAIL 0`。
- `remaining-work`: `story-E8-03` 仍未 closure-ready；Layer 3 只做到 partial coverage，`castle` / 其它内置数据集、Layer 4 Monte Carlo 和剩余测试矩阵仍待收口。

## 本轮复核证据

### 1. source-level 结构检查

执行：

```r
e <- new.env()
sys.source("lwdid-r/R/sensitivity.R", envir = e)
names(formals(e$lwdid_sensitivity))
```

结果仍显式包含：

- `max_anticipation`
- `detection_threshold`

并且对以下模式的函数体匹配均返回 `FALSE`：

- `lwdid_sensitivity <- function`
- `print.lwdid_sensitivity_comprehensive <- function`
- `summary.lwdid_sensitivity_comprehensive <- function`
- `plot.lwdid_sensitivity_comprehensive <- function`

### 2. real-data parity 回归

执行：

```bash
Rscript -e 'devtools::test(pkg="lwdid-r", filter="sensitivity-comprehensive-realdata-parity")'
```

结果：

- `PASS 14`
- `FAIL 0`
- `WARN 25`

warning 主要来自 `smoking` 数据集 `Only 1 treated unit` 的诊断提示；当前 exact / numeric parity 不以该 warning 文本一致性为 closure 条件。

### 3. oracle 状态

`_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-comparator.json`
当前记录：

- `exact_status = passed`
- `numeric_status = passed`
- `blockers = []`

## 状态建议

1. `automation-state.yaml` / `current-program.md` 不应再把结构清理写成当前主目标。
2. `parity-qa-status.yaml` 应继续把 Layer 3 表述为 partial progress，而非结构 blocker。
3. 后续优先级应回到：
   - `castle` / staggered real-data comparator；
   - Layer 4 Monte Carlo；
   - story 剩余测试矩阵。
