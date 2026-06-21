# story-E8-03 controls comparator blocker 收敛审计

## 结论

- `confirmed`: `castle` staggered / controls 路径已经不再是“待审计”状态。
  仓内现有
  `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-castle-controls-comparator.json`
  与
  `lwdid-r/tests/testthat/test-sensitivity-comprehensive-staggered-controls-parity.R`
  共同冻结了 paper-backed、bug-aware 的 controls oracle，而且
  `devtools::test(filter="sensitivity-comprehensive-staggered-controls-parity")`
  fresh rerun 当前通过。
- `confirmed`: 该 oracle 的真值不是“Python exact parity”，而是
  “论文 `(2.18)` + R simple-controls fallback + Python bug ledger waiver”。
  也就是说，`castle` controls 已经有可执行 comparator，但它是一条
  bug-aware comparator，不应再与 generic “controls comparator 未收口”混写。
- `confirmed`: controls 分支的剩余未落地 blocker 已收窄为
  `smoking` common-timing 路径的 shared frozen `X_i` fixture / comparator。
  在这条 fixture 规则写清前，`smoking` raw `X_it` controls 仍只能作为
  blocker evidence，不能写进 exact / numeric parity oracle。

## 真值锚点

### 1. `castle` controls comparator 已经具备 source-backed oracle

- qa-parity JSON
  `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-castle-controls-comparator.json`
  明确把状态写成：
  - `comparison.status = drift-confirmed`
  - `comparison.paper_backed_r_status = passed`
  - `python_reference.bug_ledger_id = PY-EST-002`
- 对应说明文档
  `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-castle-controls-comparator.md`
  明确区分：
  - R 端 paper-backed simple-controls path；
  - Python `bugged-no-controls-fallback` path；
  - no-controls reference path。

裁决：这已经不是“有没有 controls oracle”的问题，而是“oracle 是否正确按 bug
ledger 解释”的问题；当前仓内答案是肯定的。

### 2. fresh rerun 已把该 oracle 接入可执行回归

fresh 执行：

```sh
Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="sensitivity-comprehensive-staggered-controls-parity", reporter="summary")'
```

结果：`DONE`，目标套件通过。

测试文件
`lwdid-r/tests/testthat/test-sensitivity-comprehensive-staggered-controls-parity.R`
会显式检查：

- oracle 文件存在；
- `comparison.status == "drift-confirmed"`；
- `paper_backed_r_status == "passed"`；
- R warnings 中包含
  `degraded to simple controls (lw2026 eq. 2.18)`；
- 当前 controls 结果显著偏离 no-controls oracle。

裁决：`castle` controls comparator 已经从一次性 JSON 证据升级为
testthat 回归入口。

### 3. 剩余 blocker 只剩 `smoking` shared frozen `X_i`

- `Docs/lw2025.md` / `Docs/lw2026.md` 方程 `(2.17)-(2.18)`、Stata
  centered-`X` 路径、以及 Python/R 公共 validation 都要求 time-constant
  unit-level `X_i`。
- 现有 theory note
  `_automation/test-artifacts/parity/20260323-theory-parity-e8-03-smoking-controls-fixture-contract.md`
  已确认 `smoking` 常见 raw controls 是 time-varying `X_it`。

裁决：只要 shared frozen `X_i` fixture 还没明文落地，`smoking` controls
就仍是剩余 controls parity blocker；而不是 `castle`。

## 对状态文件的含义

1. `parity-manifest.yaml` 与 `parity-qa-status.yaml` 不应继续把 blocker 写成
   笼统的 “controls comparator 审计”。
2. `traceability-matrix.md` 应明确：
   `castle` controls 已有 bug-aware comparator；
   剩余 controls 缺口是 `smoking` common-timing shared frozen `X_i`。
3. `story-E8-03/tasks.md` 的后续任务备注应同步到这一级粒度，避免 controller
   误把 `castle` controls 当作未启动项。
