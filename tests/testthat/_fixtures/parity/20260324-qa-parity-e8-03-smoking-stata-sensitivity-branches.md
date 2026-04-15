# story-E8-03 smoking sensitivity branches Stata comparator

## 结论

- `confirmed`: qa-parity 本轮为 `smoking` comprehensive 的 `pre_period` 与 `no_anticipation` 两条子路径补齐了 Stata-backed executable oracle。
- `confirmed`: 当前共审计 `19` 个 `pre_period` 子规格和 `4` 个 `no_anticipation` 子规格；所有 Stata 可执行分支的 ATT 与 SE 均落在 ATT `< 1e-06`、SE `< 0.0001` 容差内。
- `bounded`: 本证据只覆盖 Task `8.4` 在 common-timing `pre_period/no_anticipation` branches 的底层 `lwdid()` 精度，不覆盖 `8.5` 的 `ri=TRUE` / RI p-value 公式漂移。
- `resolved-with-manual-oracle`: `pre_period n_pre=1` 在 Stata 端返回 `r(2001)`，但 QA contract 现已接入 theory-parity 的 paper-backed manual oracle；当前 R 的 `ATT/SE/p/CI` 与该 oracle 全部对齐，因此该单点不再阻塞 Task `8.4`。

## 数据与命令

- 数据集：`smoking`
- R 入口：`lwdid_sensitivity(type="all")`
- Stata 入口：对 `smoking.dta` 的每个 `pre_period` / `exclude_periods` 子样本分别调用 `lwdid`

## 关键覆盖

- `pre_period`：范围 `1:19`，baseline `n_pre = 19`
- `no_anticipation`：`max_anticipation_tested = 3`，`anticipation_detected = false`
- warnings：`Only 1 treated unit (N_treated=1). Results may be unreliable.`

## 裁决

- `status = passed-via-manual-oracle`
- `numeric_status = passed`
- `manual_oracle_status = passed`
- `story_task_status = closure-ready-via-manual-oracle`
- `failed_findings = 0`
- `blocked_specs = [{"branch": "pre_period", "index": 1, "stata_rc": 2001}]`
- `boundary_decisions = [{"branch": "pre_period", "index": 1, "stata_rc": 2001, "decision": "waive", "reason": "stata-implementation-boundary", "allows_story_task_closure": true, "truth_basis": ["papers-and-math", "r-python-agreement"], "evidence": ["/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260324-theory-parity-e8-03-pre-period-npre1-manual-oracle.json", "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260324-theory-parity-e8-03-pre-period-npre1-manual-oracle.md", "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260324-theory-parity-e8-03-pre-period-npre1-stata-boundary.md"]}]`

## `n_pre=1` manual oracle

- source json: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260324-theory-parity-e8-03-pre-period-npre1-manual-oracle.json`
- source status: `manual-oracle-ready`
- oracle tuple: `ATT = -0.306762820273115`, `SE = 0.084015006052754`, `p = 0.000801909718369`
- manual_oracle_findings = [{"branch": "pre_period", "index": 1, "field": "att", "status": "pass", "expected": -0.306762820273115, "actual": -0.30676282027311524, "abs_diff": 2.220446049250313e-16, "tolerance": 1e-06}, {"branch": "pre_period", "index": 1, "field": "se", "status": "pass", "expected": 0.084015006052754, "actual": 0.08401500605275375, "abs_diff": 2.498001805406602e-16, "tolerance": 0.0001}, {"branch": "pre_period", "index": 1, "field": "pvalue", "status": "pass", "expected": 0.000801909718369, "actual": 0.0008019097183691291, "abs_diff": 1.2912847874302358e-16, "tolerance": 1e-06}, {"branch": "pre_period", "index": 1, "field": "ci_lower", "status": "pass", "expected": -0.476993392318549, "actual": -0.47699339231854987, "abs_diff": 8.881784197001252e-16, "tolerance": 0.0001}, {"branch": "pre_period", "index": 1, "field": "ci_upper", "status": "pass", "expected": -0.13653224822768, "actual": -0.13653224822768065, "abs_diff": 6.38378239159465e-16, "tolerance": 0.0001}]

## 产出文件

- JSON oracle: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json`
- 当前说明: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.md`
