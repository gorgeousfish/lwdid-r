# story-E8-03 shared-data comparator evidence

## Conclusion

- `confirmed`: the shared-data Python -> R comparator for `story-E8-03` passes on the current worktree state for the synthetic common-timing fixture.
- `confirmed`: `devtools::test(filter = "sensitivity-comprehensive")` is green on the current worktree state.
- `carry-forward concern`: `lwdid-r/R/sensitivity.R` still contains repeated top-level `lwdid_sensitivity()` / comprehensive S3 method definitions, so the file remains structurally noisy even though the active wrapper now behaves correctly on this fixture.

## Evidence

1. Shared-data comparator artifact
   - report: `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-common-timing-comparator.json`
   - fixture: 40 units, 10 periods, treatment start 2006, seed 812, one CSV shared across Python and R within the comparator run
2. Python reference (`lwdid-py_v0.2.3/src/lwdid/sensitivity.py`)
   - `pre_period.sensitivity_ratio = 0.02648200276660318`
   - `transformation.difference = 0.028811007280297174`
   - `estimator.range = 0.010243749122210932`
   - no warnings
3. R current behavior (`devtools::load_all("lwdid-r")` + `lwdid_sensitivity(..., type = "all", max_anticipation = 2L, detection_threshold = 0.90)`)
   - explicit formals include `max_anticipation` and `detection_threshold`
   - returns `lwdid_sensitivity_comprehensive`
   - numeric checks against Python all pass within the comparator tolerances (`1e-6` for ATT/range/SR, `1e-4` for SE)
4. Package test signal
   - `Rscript -e "devtools::test(filter='sensitivity-comprehensive')"` returned 10 passing tests, 0 failures

## Comparator scope

- This clears Layer 2 for one shared synthetic common-timing fixture.
- This does not clear Layer 3 real-data regression, Layer 4 Monte Carlo, or the broader cleanup debt in `sensitivity.R`.

## Worktree note

- During this run, `lwdid-r/R/sensitivity.R` had file timestamp `2026-03-23 18:02:31 CST`.
- Earlier transient blocker observations from this same run were superseded by the newer worktree state; the JSON report above is the source of truth for the current state.

## Follow-up

1. Keep this comparator as the reusable Layer 2 harness for `story-E8-03`.
2. Extend the same comparator pattern to real datasets (`smoking` / `castle`) for Layer 3.
3. Re-audit whether the remaining repeated definitions in `lwdid-r/R/sensitivity.R` are harmless duplication or future drift risk.
