# story-E8-03 `smoking` frozen-controls comparator 复核

## 结论

- `confirmed`: `smoking` common-timing shared frozen `X_i` fixture 与 bug-aware comparator 已落地，不再属于“缺文件” blocker。
- `confirmed`: fresh rerun
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="sensitivity-comprehensive-common-controls-parity", reporter="summary")'`
  当前通过。
- `confirmed`: 该 comparator 的真值不是 Python exact parity，而是
  `Docs/lw2025.md` / `Docs/lw2026.md` 方程 `(2.18)`、shared frozen `X_i`
  fixture 与 `PY-EST-002` waiver 的组合。
- `confirmed`: Python 在合法 frozen `X_i` 下仍回退到
  `bugged-no-controls-fallback`；因此 `smoking` common-timing 与 `castle`
  staggered 都落在 `PY-EST-002` 的作用域内。

## 关键证据

- fixture：
  `_automation/test-artifacts/parity/e8_03_smoking_frozen_controls_fixture.csv`
- comparator：
  `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-controls-comparator.json`
- fixture note：
  `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-frozen-controls-fixture.md`
- test entry：
  `lwdid-r/tests/testthat/test-sensitivity-comprehensive-common-controls-parity.R`

## 数值与 drift 读法

- R controls demean ATT：`-0.31957385285073947`
- Python controls demean ATT：`-0.4221746008541101`
- no-controls demean ATT：`-0.42217460085410935`

裁决：Python controls demean ATT 与 no-controls 基线重合，而 R controls demean ATT
与二者存在稳定偏离；这正是 `PY-EST-002` 在 common-timing / frozen-`X_i`
路径上的 source-backed manifestation。

## 对 story 状态的含义

1. controls comparator 已从 blocker 变成 executable evidence。
2. 后续 common-controls parity 不应再追求 Python exact / numeric 一致，而应继续引用 `PY-EST-002`。
3. `story-E8-03` 的剩余未收口项回到 Task 6.3-6.5 / 8.1 / 8.2 / 8.4 / 8.5 / 8.6 与 Stata/FATAL 证据补齐。
