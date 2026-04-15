# story-E8-03 RI p-value formula source audit

## 时间

- `2026-03-24 03:26 CST`

## 目标

- 用论文、Stata、Python 与 R 源码核定 `ri=TRUE` 路径下的 RI p-value 定义是否一致。
- 判断 `story-E8-03` 当前剩余的 `FATAL-003` gap 是否仍可简单表述为
  “只差一条 story-level regression”。

## 真值锚点

1. `Docs/lw2025.md:100` 与 `Docs/lw2026.md:106` 都把 two-sided randomization
   p-value 明确写成 `c/N`，其中 `c` 是至少与观测统计量一样极端的 permutation
   次数，`N` 是总 permutation 数。
2. `lwdid_stata/lwdid.ado:436-438` 使用
   `mean(abs_res :>= abs(__b0))`，即直接计算 `count / N`。
3. `lwdid-py_v0.2.3/src/lwdid/randomization.py:477-478` 与
   `lwdid-py_v0.2.3/src/lwdid/staggered/randomization.py:527-528`
   都使用简单比例 `mean(...)`，与论文 / Stata 一致。
4. `lwdid-r/R/randomization.R:199-203` 与 `lwdid-r/R/randomization.R:482-485`
   当前都显式使用 continuity correction：
   `(1 + count) / (n_valid + 1)`。

## 裁决

- `confirmed`：R 当前 common-timing 与 staggered 两条 RI p-value 公式都偏离
  论文、Python 与 Stata 的 `c/N` 定义。
- `confirmed`：这不是 Python bug；Python 与论文 / Stata 在该点上一致，因此
  `Docs/Python包bug列表.md` 本轮无需新增条目。
- `confirmed`：该 drift 不推翻 `FATAL-003` 的架构继承结论。R 与 Stata 都已在
  RI 循环中重建 `d_perm * (X - \bar{X}_1)` 交互项，wrapper 透传 `ri=TRUE`
  的 source audit 也已经完成。
- `confirmed`：对 `story-E8-03` 而言，Task `8.5` 当前不应再被表述成
  “只差一条 `ri=TRUE` story-level regression”。更准确的说法是：
  1. 仍缺 wrapper-level `ri=TRUE` executable regression；
  2. 若后续要声称 exact RI parity，则 R 侧 `ri_pvalue` 公式本身还需要对齐论文真值。

## 对后续 parity 的操作含义

1. 在 R 尚未把 RI p-value 公式改回 `c/N` 之前，不应把当前
   `ri_pvalue` exact match 当作 source-backed oracle。
2. 后续 `qa-parity` 若先补 `FATAL-003` 的 story-level regression，建议优先比较：
   `obs_att`、`ri_distribution`、valid/failed replication 计数、warning path，
   以及 wrapper 与 direct `lwdid()` 的继承一致性。
3. 若后续 `story-worker` 选择同步修正 R 侧 RI 公式，则 common-timing 与
   staggered 两条实现应一起对齐，避免 Layer 5 再次出现内部双轨定义。
