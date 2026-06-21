# story-E8-03 `pre_period n_pre=1` Stata boundary source audit

## 时间

- `2026-03-24 04:26 CST`

## 要回答的问题

- `smoking` comprehensive `pre_period n_pre=1` 在 R 端可执行、Stata 端返回 `r(2001)`；
  该现象是论文定义不允许、R/Python 偏离论文，还是 Stata ado 的实现边界？

## 读取的 source

1. `Docs/lw2025.md:42` 明确 common-timing 识别只要求“至少一个” pre-treatment period。
2. `Docs/lw2025.md:180-184` 与 `Docs/lw2026.md:186-194` 明确给出单期 reference
   transformation `Y_it - Y_i,S-1`（公式 `(2.21)`），并说明可以只使用 treatment
   前最后一期，乃至使用任意 pre-treatment average。
3. `lwdid-r/R/sensitivity.R:223-226` 把 `demean` / `demeanq` 的
   `min_pre` 固定为 `1L`。
4. `lwdid-py_v0.2.3/src/lwdid/sensitivity.py:983-989` 同样把 `demean` /
   `demeanq` 的 `min_pre` 固定为 `1`。
5. `lwdid_stata/lwdid.ado:243-246` 在 `rolling(demean)` 下对每个 unit 直接运行
   `regress y if id==ii & post==0`；当某 unit 仅剩 1 条 pre-period 观测时，
   Stata `regress` 会直接报 `insufficient observations`。

## 复现

### Stata 端最小复现

- 数据：`/Users/cxy/Desktop/lwdid_r/lwdid_stata/smoking.dta`
- 过滤：`keep if post == 1 | year == 1988`
- 命令：`lwdid lcigsale d, ivar(state) tvar(year) post(post) rolling(demean)`
- 观测到的失败文本：
  - `insufficient observations`
  - `RC=2001`

### R / Python 对照

- R `lwdid_sensitivity_pre_period()` 与 Python `robustness_pre_periods()` 都会把
  common-timing `demean` 的 pre-period search range 从 `1` 开始。
- qa-parity 现有 oracle
  `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json`
  已确认 `n_pre=2:19` 与 Stata 一致；`n_pre=1` 是唯一不可执行分支。

## 裁决

- `confirmed`: 论文真值允许 common-timing `demean` 使用单个 pre-treatment period。
  `lw2025/lw2026` 的 `(2.21)` 直接对应 `Y_it - Y_i,S-1`，因此 `n_pre=1`
  不是 theory violation。
- `confirmed`: 当前 R 与 Python 在 `n_pre=1` 上没有对论文的 drift；两边的
  auto-detected minimum pre-period 都是 `1`。
- `confirmed`: `smoking` 的 `pre_period n_pre=1` 失败来自
  `lwdid_stata/lwdid.ado` 的 unit-wise `regress` 实现路径，而不是论文定义、
  也不是 R/Python 的实现错误。该边界应归类为 `Stata implementation boundary`。

## 对 story-E8-03 的含义

- Task `8.4` 不应继续被表述为“R 侧仍有 `pre_period n_pre=1` theory / numeric drift”。
- 若 closure 要求所有 `pre_period` branches 都有 Stata-backed oracle，则
  `n_pre=1` 需要二选一：
  1. 手工按 `(2.21)` / cross-sectional ATT 公式重建 oracle；
  2. 明确登记为 Stata ado 无法执行的 waiver/boundary。
- 与此同时，`_automation/test-artifacts/parity/20260324-controller-e8-03-namespace-dispatch-drift.md`
  记录的 namespace / S3 dispatch drift 仍是独立 blocker，不受本裁决影响。
