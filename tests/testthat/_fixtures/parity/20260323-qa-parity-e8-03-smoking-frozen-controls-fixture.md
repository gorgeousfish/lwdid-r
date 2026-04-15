# story-E8-03 `smoking` frozen-controls fixture 证据

## 结论

- `decision`: 对 `smoking` common-timing controls comparator 采用 shared frozen `X_i` fixture。
- `rule`: 对每个 `state`、每个 control 列，仅使用 `post == 0` 的观测做 state-level mean。
- `why`: 该规则显式避免 post-treatment leakage，并把 raw `X_it` 冻结为 comparator 可用的 unit-level `X_i`。
- `scope`: 这是 parity fixture 规则，不是把 R 内部 pre-period proxy helper 直接提升为公共 API 真值。

## 生成窗口

- 起始年份：`1970`
- 结束年份：`1988`
- 州数量：`39`

## 非缺失覆盖

- `lnincome`: 每州 pre-period 非缺失条数固定为 `17`
- `beer`: 每州 pre-period 非缺失条数固定为 `5`
- `age15to24`: 每州 pre-period 非缺失条数固定为 `19`
- `lretprice`: 每州 pre-period 非缺失条数固定为 `19`

## 产物

- fixture csv: `/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e8_03_smoking_frozen_controls_fixture.csv`

## 约束

- comparator 与 testthat 回归都必须直接复用该 fixture，不能让 Python 或 R 各自临时再聚合一遍。
- raw `smoking` controls 仍然不是合法公共 comparator 输入；共享 fixture 只是为了冻结 bug-aware 审计边界。
