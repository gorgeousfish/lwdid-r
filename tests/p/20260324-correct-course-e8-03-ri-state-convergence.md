# story-E8-03 RI / state convergence audit

## 时间

- `2026-03-24 03:58 CST`

## fresh rerun

- 命令：
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="randomization|sensitivity-comprehensive|sensitivity-comprehensive-fatal-inheritance-parity|sensitivity-comprehensive-stata-precision-parity", reporter="summary")'`
- 结果：目标 suites 当前全部通过；仅保留既有 warning，
  未出现新的 error / failure。

## 复核到的代码事实

1. `lwdid-r/R/randomization.R` 的 common / staggered 两条 RI p-value
   公式当前都已改为 simple proportion `c / N`，不再使用 continuity correction。
2. `lwdid-r/tests/testthat/test-randomization.R` 已同步把 formula oracle
   对齐到 `c / N`。
3. `lwdid-r/tests/testthat/test-sensitivity-comprehensive-fatal-inheritance-parity.R`
   当前已包含 common / staggered `type="transformation"` + controls +
   `ri=TRUE` 的 story-level regression，并显式断言 wrapper 子调用返回的
   `ri_pvalue` 满足 `c / N`。
4. `20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json`
   继续表明 Task `8.4` 当前仅剩 `pre_period n_pre=1` 的 Stata `r(2001)` 边界；
   `pre_period n_pre=2:19` 与 `no_anticipation exclude=0:3` 均已有
   executable Stata-backed ATT/SE 锚点。

## 裁决

- `confirmed`: Task `8.5` 已从 state / parity blocker 转为 executable evidence。
- `confirmed`: `automation-state.yaml`、`current-program.md`、
  `parity-manifest.yaml`、`parity-qa-status.yaml` 与
  `traceability-matrix.md` 不应继续把 `8.5` 表述为当前 blocker。
- `confirmed`: 当前 active story 的直接缺口已收敛为：
  - Task `8.4` 中仅剩的 `pre_period n_pre=1` Stata `r(2001)` 边界；
  - `E8-04` / `E8-05` 所需 `R/clustering_diagnostics.R`、
    `R/selection_diagnostics.R` 及对应 testthat 入口仍不存在。
- `not-needed`: 本轮无需更新 `Docs/Python包bug列表.md`，因为本次 drift 属于
  R-side state convergence，而非新的 Python bug。
