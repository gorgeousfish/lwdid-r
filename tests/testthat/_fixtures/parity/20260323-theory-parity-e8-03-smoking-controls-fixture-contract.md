# story-E8-03 `smoking` common-timing controls comparator 输入契约审计

## 结论

- `confirmed`: `Docs/lw2026.md` 方程 `(2.17)` 与 `(2.18)`、`Docs/lw2025.md`
  对应 conditional PT 表述，以及 `lwdid_stata/lwdid.ado` 的 centered-`X`
  路径都把 controls 定义为 time-constant unit-level `X_i`，不是 raw panel
  `X_it`。
- `confirmed`: `smoking.csv` 的常见 raw controls `lnincome`、`beer`、
  `age15to24`、`lretprice` 在 39 个州内都随时间变化；full-sample
  within-state `max sd` 分别为 `0.1756441`、`6.3332412`、`0.0196662`、
  `0.6825728`。
- `confirmed`: Python 与 R 的公共 common-timing validation 都会拒绝把这些
  raw columns 当作合法 `controls` 传入 `lwdid_sensitivity()`；因此当前没有
  source-backed 的 `smoking` controls comparator 可直接执行。
- `boundary`: `lwdid-r/R/lwdid.R` 内部 `.extract_controls()` 虽保留
  “pre-period mean proxy” 逻辑，但公共入口会先经过
  `validate.R::.validate_time_invariant_controls()`；该 proxy 不是当前公共 API
  的 parity 真值，也不是论文显式规定的 frozen-`X_i` 规则。
- `decision`: 后续若要补 `smoking` controls parity，必须先产出 shared frozen
  `X_i` fixture，并把生成规则写成单独证据；在此之前，`smoking` raw controls
  只能作为 blocker evidence，不能写入 comparator oracle。

## 真值锚点

### 1. 论文中的 controls 对象是 unit-level `X_i`

- `Docs/lw2026.md:154-164`：
  conditional PT 写成
  `E[Y_it(0)-Y_i1(0) | D_i, X_i] = E[Y_it(0)-Y_i1(0) | X_i]`，并在
  `(2.18)` 中使用
  `Delta Ybar_i = alpha + tau D_i + X_i beta + U_i`。
- `Docs/lw2025.md` 对应 conditional PT / transformed-regression 表述同样以
  `X_i` 为对象，而不是逐期 `X_it`。

裁决：只要要做 paper-backed common-timing controls comparator，输入对象就必须是
每个 unit 一行的固定 covariate。

### 2. Stata 参考实现也假定 `X` time-invariant

- `lwdid_stata/lwdid.ado:320` 明确写着
  `build centered X and interactions ONCE (time-invariant X)`。
- `lwdid_stata/lwdid.ado:346` 进一步说明
  `Since X is time-invariant (X_it = X_i), we need to average over units, not observations`。

裁决：Stata 侧并没有提供“直接接收 raw `X_it` 再自动冻结”的 comparator 真值。

## 公共实现证据

### 3. Python 公共 validation 会拒绝 `smoking` raw controls

fresh 复跑：

- 调用：
  `sensitivity_analysis(..., controls=['lnincome','beer','age15to24','lretprice'])`
- 结果：
  comprehensive wrapper 继续返回对象，但各 specification / exclusion /
  transformation path 都发出同一类 warning：
  `Control variables must be time-invariant ... requires time-constant controls X_i, not time-varying X_it`
- Python 同时给出聚合建议：
  first period / unit mean / pre-treatment value / domain-appropriate method。

裁决：Python 并未把 raw `X_it` 视为合法 comparator 输入；它只是提示未来需要先做
冻结。

### 4. R 公共 validation 与 Python 一致

fresh 复跑：

- 调用：
  `lwdid_sensitivity(..., controls = c('lnincome','beer','age15to24','lretprice'))`
- 结果：
  wrapper 返回 comprehensive 对象，但 pre-period / transformation 等子路径均报
  `Control variables must be time-invariant (constant within each unit)`，
  并列出 39 个州内的 within-state 波动。
- 源码位置：
  `lwdid-r/R/validate.R:678-718` 的
  `.validate_time_invariant_controls()` 在 common-timing mode 下强制拦截时变
  controls；
  `lwdid-r/R/validate.R:2142-2156` 显示该校验发生在估计入口之前。

裁决：R 的公共真值也不是“接受 raw `X_it` 然后静默聚合”。

## 边界说明：`.extract_controls()` 不能替代 fixture 规范

- `lwdid-r/R/lwdid.R:636-683` 的 `.extract_controls()` 会对时变 controls 发出
  `pre-period means used as time-invariant proxy` warning，并在 pre window 上求均值。
- 但 `lwdid-r/R/validate.R:2142-2156` 先执行
  `.validate_time_invariant_controls()`，所以公共调用链在 common-timing mode 下
  不会把 raw `smoking` controls 放进 `.extract_controls()` 后再继续估计。

裁决：

- 这个 proxy 更像内部降级草稿，不是当前对外 API 的 source-backed parity 规则。
- 后续如要采用 “pre-period mean freeze” 作为 comparator fixture 方案，必须由新证据
  明确宣布，而不能把该 helper 的存在直接当作理论或接口承诺。

## 对后续 comparator 的约束

1. `smoking` controls comparator 必须先定义 frozen `X_i` 生成规则。
2. 该规则必须在 Python / R 两侧共用同一 fixture，不能让任一端自行聚合。
3. 在未形成 shared frozen fixture 前，不得把 `smoking` raw controls 写进
   `story-E8-03` 的 exact / numeric parity oracle。
4. `castle` controls drift 与本条不同：`castle` 的 controls 是合法 `X_i`，
   它的 blocker 来自 Python 小样本 `eq. 2.18` fallback 漂移，已由 bug ledger
   单独覆盖。
