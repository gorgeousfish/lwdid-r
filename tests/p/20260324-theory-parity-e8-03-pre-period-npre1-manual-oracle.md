# story-E8-03 `pre_period n_pre=1` manual oracle

## 时间

- `2026-03-24 05:32 CST`

## Target

- 为 `smoking` comprehensive `pre_period n_pre=1` 生成不依赖 Stata `lwdid.ado` 的 paper-backed ATT/SE oracle。
- 判断该单点当前是否仍属于“缺少数值真值”，还是已经具备可供 qa-parity / story-worker 接入的 manual oracle。

## Invariant Object

- 共同时点、`rolling(demean)`、仅使用最后一个 pre-treatment period 的 transformed outcome。
- 论文锚点：
  - `Docs/lw2025.md:180-182`
  - `Docs/lw2026.md:186-188`
- 对应公式：
  - `(2.21)`：`\mathring{Y}_{it} = Y_{it} - Y_{i,S-1}`
  - `(2.20)`：对 transformed outcome 做横截面 ATT 回归
  - `(2.10)`：小样本 `t` 推断

## Assumptions And Notation

- 数据：`smoking`
- design：common timing
- reference period：`S-1 = 1988`
- post periods：`1989:2000`
- treated units：`N_1 = 1`
- control units：`N_0 = 38`
- total units：`N = 39`
- 无 controls；与现有 `pre_period` branch 的 `n_pre=1` 规格一致

## Main Derivation

1. 由论文 `(2.21)`，对每个单位、每个处理后年份定义
   `\mathring{Y}_{it} = Y_{it} - Y_{i,1988}`。
2. 对每个单位取 post-period average：
   `\overline{\mathring{Y}}_i = (1/12) * \sum_{t=1989}^{2000} (Y_{it} - Y_{i,1988})`。
3. 对横截面样本 `{(\overline{\mathring{Y}}_i, D_i)}` 运行 `\overline{\mathring{Y}}_i on 1, D_i`。
4. `D_i` 的 OLS 系数就是 `n_pre=1` branch 的 paper-backed ATT；标准误和区间按同一回归获得。

## Oracle Values

- `ATT = -0.306762820273115`
- `SE = 0.084015006052754`
- `p = 0.000801909718369`
- `95% CI = [-0.476993392318549, -0.136532248227680]`

## Cross-Implementation Check

### Manual vs current R

- R current (`pkgload::load_all()` + `lwdid_sensitivity(type="all")`) 的 `n_pre=1` 规格为：
  - `ATT = -0.306762820273115`
  - `SE = 0.084015006052754`
  - `p = 0.000801909718369`
  - `95% CI = [-0.476993392318550, -0.136532248227681]`
- 与 manual oracle 的绝对差：
  - `ATT = 0`
  - `SE = 0`
  - `p = 0`
  - `CI lower/upper ~= 9.99e-16`

### Manual vs current Python

- Python `robustness_pre_periods(..., rolling="demean")` 的 `n_pre=1` 规格为：
  - `ATT = -0.306762808333334`
  - `SE = 0.084014997113155`
  - `p = 0.000801909154138`
  - `95% CI = [-0.476993362265420, -0.136532254401247]`
- 相对 manual oracle 的绝对差：
  - `ATT = 1.193978099722770025e-08`
  - `SE = 8.939599005497456119e-09`
  - `p = 5.642310000296413275e-10`
  - `CI lower = 3.005312898540779543e-08`
  - `CI upper = 6.173566990952394917e-09`
- 以上均落在现有 parity 容差 `ATT < 1e-6`、`SE < 1e-4` 内。

## Decision

- `confirmed`: `n_pre=1` 已具备 paper-backed manual oracle，不再属于“缺少理论或数值真值”。
- `confirmed`: Stata `r(2001)` 仍然只是 `lwdid.ado` 的实现边界，不影响该 branch 的论文合法性。
- `confirmed`: 当前未发现新的 Python-paper drift；因此本轮不更新 `Docs/Python包bug列表.md`。

## Implication For Task 8.4

- Task `8.4` 的剩余动作已不再是排查 R/Python 公式。
- 下一步只剩二选一：
  1. qa-parity / story-worker 将本 manual oracle JSON 接入 `n_pre=1` branch 的 ATT/SE 断言；
  2. controller 明确保留 “Stata branch unavailable, use manual oracle / waiver” 的 closure 文案。
