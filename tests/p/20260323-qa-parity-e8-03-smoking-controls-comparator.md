# story-E8-03 `smoking` controls comparator bug-aware oracle

## 结论

- `status`: `drift-confirmed`
- `paper-backed R`: `passed`。当前 R 端在 shared frozen `X_i` fixture 上继续走论文 `eq. 2.18` simple-controls fallback。
- `Python`: `bugged-no-controls-fallback`。即使输入已冻结为合法 `X_i`，Python 仍因 `N_1 <= K+1` 直接忽略 controls。

## frozen `X_i` fixture

- fixture artifact: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e8_03_smoking_frozen_controls_fixture.csv`
- rule: `state-level pre-treatment mean over observed rows with post == 0`
- window: `1970 - 1988`

## 关键数值

- R controls demean ATT: `-0.319573852850739`
- Python controls demean ATT: `-0.422174600854110`
- no-controls demean ATT: `-0.422174600854109`
- R controls vs no-controls |demean diff|: `0.102600748003370`
- R controls rel_range: `0.321054889466474`

## 关键 warning

- R 端包含 `degraded to simple controls (lw2026 eq. 2.18)`，说明当前走的是 paper-backed fallback。
- Python 端包含 `Controls not applied ...`，说明当前仍直接丢弃 controls。

## 证据路径

- contract audit: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260323-theory-parity-e8-03-smoking-controls-fixture-contract.md`
- fixture evidence: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-frozen-controls-fixture.md`
- bug ledger: `/Users/cxy/Desktop/lwdid_r/Docs/Python包bug列表.md`
- no-controls oracle: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-comparator.json`

## R warning 摘要

- Only 1 treated unit (N_treated=1). Results may be unreliable.
- Insufficient sample for full interaction model (need N_1>5 and N_0>5, got N_1=1, N_0=38), degraded to simple controls (lw2026 eq. 2.18)
- glm.fit: algorithm did not converge
- glm.fit: fitted probabilities numerically 0 or 1 occurred
- glm did not converge. PS coefficients may be unreliable.
