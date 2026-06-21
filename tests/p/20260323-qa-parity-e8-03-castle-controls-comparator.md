# story-E8-03 `castle` controls comparator bug-aware oracle

## 结论

- `status`: `drift-confirmed`
- `paper-backed R`: `passed`。当前 R 端继续走论文 `eq. 2.18` simple-controls fallback。
- `Python`: `bugged-no-controls-fallback`。controls 路径数值回退到既有 no-controls oracle，不应当作 parity 真值。

## 关键数值

- R controls demean ATT: `0.064562798306184`
- Python controls demean ATT: `0.091745387139613`
- no-controls demean ATT: `0.091745387139613`
- R controls vs no-controls |demean diff|: `0.027182588833429`

## 关键 warning

- R 端包含 `degraded to simple controls (lw2026 eq. 2.18)`，说明当前走的是 paper-backed fallback。
- Python 端包含 `Controls not included ...`，说明当前直接丢弃了 controls。

## 证据路径

- theory audit: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260323-theory-parity-e8-03-controls-small-sample-drift.md`
- bug ledger: `/Users/cxy/Desktop/lwdid_r/Docs/Python包bug列表.md`
- no-controls oracle: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260323-qa-parity-e8-03-castle-comparator.json`

## R warning 摘要

- Unbalanced panel: obs per unit ranges from 3 to 7. 20 units (40.0%) incomplete.
Diagnostics: 100.0% usable for demean, 0.0% usable for detrend.
- aggregate='cohort' with control_group='not_yet_treated': auto-switching to control_group='never_treated' for aggregation. Cross-cohort aggregation requires a common control group (NT units).
- Insufficient sample for full interaction model (need N_1>4 and N_0>4, got N_1=1, N_0=29), degraded to simple controls (lw2026 eq. 2.18)
- Insufficient sample for full interaction model (need N_1>4 and N_0>4, got N_1=4, N_0=29), degraded to simple controls (lw2026 eq. 2.18)
- Insufficient sample for full interaction model (need N_1>4 and N_0>4, got N_1=2, N_0=29), degraded to simple controls (lw2026 eq. 2.18)
- Cohort g=2005: controls degraded (N_treat=1 or N_ctrl=29 insufficient, need >4)
- Cohort g=2007: controls degraded (N_treat=4 or N_ctrl=29 insufficient, need >4)
- Cohort g=2008: controls degraded (N_treat=2 or N_ctrl=29 insufficient, need >4)
- Cohort g=2009: controls degraded (N_treat=1 or N_ctrl=29 insufficient, need >4)
- Unbalanced panel: obs per unit ranges from 4 to 8. 20 units (40.0%) incomplete.
Diagnostics: 100.0% usable for demean, 100.0% usable for detrend.
- Unbalanced panel: obs per unit ranges from 5 to 9. 20 units (40.0%) incomplete.
Diagnostics: 100.0% usable for demean, 100.0% usable for detrend.
- Unbalanced panel: obs per unit ranges from 6 to 10. 20 units (40.0%) incomplete.
Diagnostics: 100.0% usable for demean, 100.0% usable for detrend.
- Unbalanced panel: obs per unit ranges from 7 to 11. 20 units (40.0%) incomplete.
Diagnostics: 100.0% usable for demean, 100.0% usable for detrend.
- Unbalanced panel: obs per unit ranges from 10 to 11. 21 units (42.0%) incomplete.
Diagnostics: 100.0% usable for demean, 100.0% usable for detrend.
- Unbalanced panel: obs per unit ranges from 9 to 11. 21 units (42.0%) incomplete.
Diagnostics: 100.0% usable for demean, 100.0% usable for detrend.
- Unbalanced panel: obs per unit ranges from 8 to 11. 21 units (42.0%) incomplete.
Diagnostics: 100.0% usable for demean, 100.0% usable for detrend.
- IPW weight CV = 2.13 > 2.0. Possible overlap violation.
- glm.fit: algorithm did not converge
- glm.fit: fitted probabilities numerically 0 or 1 occurred
- glm did not converge. PS coefficients may be unreliable.
- no non-missing arguments to max; returning -Inf
- NAs introduced by coercion to integer range
- Small treated sample (4). IPWRA estimates may be unstable.
- IPWRA weight CV = 2.13 > 2.0. Possible overlap violation.
- 48.5% of PS < 0.05 and 0.0% > 0.95. Overlap assumption may be violated.
- Small treated sample (2). IPWRA estimates may be unstable.
- 74.2% of PS < 0.05 and 0.0% > 0.95. Overlap assumption may be violated.
- 8/20 (g,r) pairs skipped (estimation_error: 8)
