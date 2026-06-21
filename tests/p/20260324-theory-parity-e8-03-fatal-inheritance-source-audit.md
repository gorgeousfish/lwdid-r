# story-E8-03 FATAL 继承 source audit

## 时间

- `2026-03-24 02:30 CST`

## 目标

- 用论文与源代码核定 `story-E8-03` 的 comprehensive sensitivity wrapper
  是否已经继承底层 `lwdid()` 的 FATAL 防护。
- 把 Task `8.5` 从“理论是否明确”收缩为“还缺哪类可执行回归”。

## 真值锚点

1. `Docs/lw2025.md` 式 `(4.8)` 明确 staggered not-yet-treated 对照组必须满足
   `G_i > r`，不能把 `G_i = r` 的单位纳入 control pool。
2. `Docs/lw2026.md` 式 `(7.10)` 把 cohort-level aggregation 的回归样本限定为
   `D_{ig} + D_{i\infty} = 1`，因此 cohort / overall aggregation 的共同对照组
   必须是 never-treated。
3. `Docs/lw2026.md` 式 `(7.18)` 要求 never-treated 单位使用
   `\sum_g \hat{\omega}_g \overrightarrow{Y}_{ig}` 的加权平均，而不是单一 cohort
   transform。
4. `lwdid-r/R/control_groups.R` 已把 FATAL-001 写成 source comment +
   executable mask：`not_yet_treated` 只能用 strict `g_vals > r`。
5. `lwdid-r/R/lwdid.R` 已把 cohort / overall aggregation 下的 NT-only 约束写成
   FATAL-004：若当前 control group 不是 `never_treated`，则自动切换并重跑
   staggered estimation。
6. `lwdid-r/R/aggregate.R` 已把 FATAL-002 写成 NT 单位对所有成功 cohort
   transform 的加权平均，并包含缺失时的 per-unit renormalization。
7. `lwdid-r/R/randomization.R` 已把 controls 路径的 RI 设计矩阵写成
   “固定 `x_c`、逐次重建 `d_perm * x_c` interaction”，即 FATAL-003 的实现位点。

## wrapper source audit

### 1. Comprehensive wrapper 没有重写 FATAL 逻辑

- `lwdid-r/R/sensitivity.R` 的 `lwdid_sensitivity_pre_period()`、
  `sensitivity_no_anticipation()`、`.run_transformation_comparison()`、
  `.run_estimator_comparison()` 都通过 `do.call(lwdid, ...)` 执行底层估计。
- 这些入口统一使用 `.filter_sensitivity_args(extra_args, lwdid)`，
  所以 `control_group`、`aggregate`、`ri`、`rireps`、`seed` 等底层参数都会继续传给
  `lwdid()`，而不是在 wrapper 层被吞掉。

### 2. Task 8.5 的理论对象已经明确

- `FATAL-001`：
  comprehensive wrapper 若在 staggered 场景下传入 `control_group="not_yet_treated"`，
  真值仍然是 `G_i > r` 的 strict mask，而不是 `G_i >= r`。
- `FATAL-004`：
  comprehensive wrapper 若在 staggered 场景下使用默认 `aggregate="cohort"`，
  或显式传入 `aggregate="overall"`，则底层 `lwdid()` 仍会强制使用
  never-treated 共同对照组，并在需要时发出 auto-switch warning。
- `FATAL-002`：
  上述 aggregation 一旦发生，NT 单位的 aggregated outcome 仍由
  `aggregate.R` 里的 weighted-average engine 构造；wrapper 没有单独拼装
  aggregated counterfactual，因此不存在另一套 sensitivity-only 数学定义。

### 3. FATAL-003 的位置要单独说明

- Task `8.5` 当前文字只点名 `FATAL-001` 与 `FATAL-004`，但 REQ-38 写的是
  “底层 `lwdid()` 自动继承四大 FATAL 防护”。
- source audit 显示：如果 comprehensive sensitivity 调用显式传入 `ri=TRUE`，
  wrapper 会把 `ri` 继续传入 `lwdid()`，最终走到 `randomization_inference()` /
  `ri_staggered()`，因此 FATAL-003 在架构上是“已继承”的。
- 但当前 story 的已落地回归矩阵并没有针对 `ri=TRUE` 的 comprehensive 路径给出
  executable oracle，所以 FATAL-003 现在仍应记作“未做 story 级回归覆盖”，
  而不是“理论未定”。

## 裁决

- `confirmed`：Task `8.5` 当前不是公式来源不清的问题；论文、`lwdid.R`、
  `control_groups.R`、`aggregate.R`、`randomization.R` 已经把 FATAL-001/002/003/004
  的真值和实现位点写清楚。
- `confirmed`：`story-E8-03` comprehensive wrapper 已通过参数透传继承这些防护；
  剩余缺口是 executable regression，要验证 warning / numeric oracle 是否在 wrapper
  层可观察，而不是再写一轮理论判词。
- `confirmed`：Task `8.5` 的最小可执行收口应至少覆盖：
  1. staggered + `control_group="not_yet_treated"` 的 strict-mask 数值差异或 warning；
  2. staggered + `aggregate="cohort"` / `"overall"` 的 NT-only auto-switch 继承；
  3. 至少一条 aggregation oracle，证明 comprehensive wrapper 没有绕开
     `eq. (7.18)` 的 NT weighted-average 构造。
- `suspected`：若要把 REQ-38 的 “四大 FATAL 全继承” 一次性打勾，仍需追加
  `ri=TRUE` 的 story-level regression；否则 Task `8.5` 最好明确写成
  “当前收口目标先锁定 FATAL-001/002/004，FATAL-003 另列 RI 覆盖”。

## 对当前 story 的含义

1. `8.5` 继续保持 open，但应被视为 regression-gap，而不是 theory-gap。
2. 后续 `story-worker` / `qa-parity` 若继续实现 `8.5`，优先补 staggered executable
   regression，不要再回头争论 `G_i > r`、NT-only aggregation 或 `eq. (7.18)` 的含义。
3. `traceability` / `parity-manifest` 应记录：Layer 5 仍缺 FATAL inheritance
   regression，但 source-backed theoretical boundary 已收口。
