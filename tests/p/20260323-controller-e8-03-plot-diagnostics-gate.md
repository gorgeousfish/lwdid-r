# story-E8-03 implementation gate after plot diagnostics

## 结论

- `confirmed`: `story-E8-03` 仍不能离开 `implementation` phase。
- `confirmed`: `sensitivity.R` 的综合入口 / S3 结构漂移已解除，`smoking`
  real-data parity 也已通过；当前 gate blocker 不再是结构问题。
- `confirmed`: fresh `plot-diagnostics` 仍存在可重复的实现回归，因此 controller
  本轮保持 `active_story=story-E8-03` 不变，并把当前目标收敛到
  sensitivity 可视化回归 + 剩余 parity。

## 本轮验证

### 1. positive signal

执行：

```bash
Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="sensitivity-comprehensive", reporter="summary")'
```

结果：

- `sensitivity-comprehensive: PASS`
- `sensitivity-comprehensive-realdata-parity: PASS`

这说明 comprehensive 入口与 `smoking` oracle 回归当前可运行。

### 2. gate blocker

执行：

```bash
Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="plot-diagnostics", reporter="summary")'
```

结果显示至少以下失败仍存在：

- TC-10.2.2：预期效应图缺 `GeomLine`
- TC-10.2.9：`show_threshold=FALSE` 未移除 threshold band
- TC-10.2.23：规格曲线 threshold band 仍画成 `1.5-2.5`，而非
  基准值 `2.0` 的 `±25%` 即 `1.8-2.2`
- TC-10.2.24 / TC-10.2.25 / TC-10.2.33：`show_ci=FALSE` 时 errorbar
  仍保留，且 comprehensive 未正确向子图透传
- TC-10.2.27：预期效应图缺推荐排除点的 `GeomVline`
- TC-10.2.31 / TC-10.2.32：unconverged 路径仍直接抛 error，而不是
  warning + empty ggplot

## 状态建议

1. `phase` 保持 `implementation`。
2. `active_story` 保持 `story-E8-03`。
3. 先由 `story-worker` 修 sensitivity plot regressions，再由 `qa-parity`
   扩展 `castle` / Layer 4 Monte Carlo。
