# story-E8-03 controls blocker 状态收敛复核

## 时间

- `2026-03-24 00:28 CST`

## 目标

- 校验 `story-E8-03` 的 common-timing controls 路径是否仍把
  `smoking` frozen-controls comparator 当作当前 blocker。
- 以论文与数学原理为真值，确认当前应追踪的对象是
  `Docs/lw2026.md` 方程 `(2.18)` 下的 simple-controls ATT，而不是 Python
  小样本 `bugged-no-controls-fallback`。

## 真值锚点

1. `Docs/lw2026.md` 方程 `(2.17)` 与 `(2.18)` 把 controls 定义为
   time-invariant unit-level `X_i`，并允许在 `N > K + 2` 时估计
   `\Delta \bar{Y}_i = \alpha + \tau D_i + X_i \beta + U_i`。
2. `lwdid_stata/lwdid.ado` 的 centered-`X` 路径也假定 `X_it = X_i`。
3. `Docs/Python包bug列表.md` 的 `PY-EST-002` 已确认 Python 在合法 small-sample
   controls 下仍会错误丢弃 paper-valid simple-controls 回归。

## fresh verification

- 回归命令：
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="sensitivity-comprehensive-common-controls-parity", reporter="summary")'`
- 结果：`passed`
- 直接证据：
  - comparator:
    `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-controls-comparator.json`
  - shared frozen `X_i` fixture:
    `_automation/test-artifacts/parity/e8_03_smoking_frozen_controls_fixture.csv`
  - test entry:
    `lwdid-r/tests/testthat/test-sensitivity-comprehensive-common-controls-parity.R`
  - waiver note:
    `_automation/test-artifacts/parity/20260323-theory-parity-e8-03-smoking-controls-waiver-resolution.md`

## 裁决

- `confirmed`: `smoking` common-timing frozen-controls comparator 当前仍然可执行，
  不再是缺文件或缺 oracle blocker。
- `confirmed`: 当前 parity 真值是 “shared frozen `X_i` fixture + paper-backed
  `eq. 2.18` simple-controls + `PY-EST-002` waiver”，不是 Python exact parity。
- `confirmed`: `automation-state.yaml` 与 `current-program.md` 若仍把
  `20260323-qa-parity-e8-03-smoking-controls-comparator.json` /
  `e8_03_smoking_frozen_controls_fixture.csv` 记成未落地，则属于状态漂移。

## 对当前 story 的含义

1. `story-E8-03` 的剩余缺口不再包含 generic controls comparator。
2. 当前应继续追踪的缺口回到 Task `6.3-6.5 / 8.1 / 8.2 / 8.4 / 8.5 / 8.6`。
3. `E8-04` / `E8-05` 所需的
   `lwdid-r/R/clustering_diagnostics.R`、
   `lwdid-r/R/selection_diagnostics.R` 及对应 testthat 入口仍不存在，因此 active
   story 继续保持 `story-E8-03`。
