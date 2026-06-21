# story-E8-03 plot diagnostics source-backed oracle audit

## 结论

- `confirmed`: `TC-10.2.23` 当前把规格曲线 threshold band 解释为
  `baseline_att ± robustness_threshold * |baseline_att|`，但
  `story-E8-03` 的 requirements/design 与 Python 参考实现都把该图层固定为
  `±25% of baseline ATT`。因此这条失败目前属于 test oracle drift，而不是
  source-backed parity failure。
- `confirmed`: `TC-10.2.27` 当前要求 no-anticipation 图包含 `GeomVline`，
  但 requirements/design 与 Python 参考实现都只要求在推荐排除点绘制红色空心圆。
  因此“缺少竖线”不应被视为论文/Python parity blocker。
- `confirmed`: `TC-10.2.2` 对 no-anticipation 图包含折线层
  (`GeomLine`) 的要求与 Python `fmt='o-'` 一致，仍属于 source-backed
  实现回归。

## 真值锚点

### 1. 论文只提供稳健性分析动机，不提供新增图形 oracle

- `Docs/lw2026.md:627-631` 只说明可通过改变 pre-treatment periods
  研究估计稳健性。
- `Docs/lw2025.md:451-455` 只说明若怀疑 no anticipation 失效，可排除处理前若干期
  重新估计。

结论：plot 级别的具体元素必须由 story spec 和参考实现共同锚定，不能从当前 failing
 tests 自行反推。

### 2. story spec 对 plot contract 的明示要求

- `.kiro/specs/story-E8-03/requirements.md:90-94`
  - `REQ-29`: 规格曲线图包含 `±25%` 稳健性带；
  - `REQ-31`: 预期效应图在推荐排除点使用红色空心圆标记。
- `.kiro/specs/story-E8-03/design.md:336-349`
  明确 no-anticipation 图的推荐标记为红色空心大圆。
- `.kiro/specs/story-E8-03/design.md:318-322`
  明确稳健性带为 `baseline ATT ± 25% * |baseline ATT|`。

### 3. Python 参考实现

- `lwdid-py_v0.2.3/src/lwdid/sensitivity.py:479-519`
  的 pre-period plot 在 baseline 附近使用 `0.25 * abs(self.baseline_spec.att)`
  构造 robustness band。
- `lwdid-py_v0.2.3/src/lwdid/sensitivity.py:676-713`
  的 no-anticipation plot 用 `fmt='o-'` 绘制点线，并在推荐排除点使用
  `ax.scatter(..., facecolors='none', edgecolors='red')`；没有 `axvline`
  或其它竖线标记。

## 当前 R 与 tests 的对应关系

### R 当前实现

- `lwdid-r/R/sensitivity.R:857-866`
  在 no-anticipation 图里只画 errorbar、point 和推荐点空心圆，确实缺少 Python
  对标的折线层。
- `lwdid-r/R/sensitivity.R:896-908`
  以固定 `0.25 * abs(baseline_att)` 构造规格曲线 threshold band，
  这与 requirements/design/Python 一致。

### 当前 test oracle 漂移

- `lwdid-r/tests/testthat/test-plot-diagnostics.R:292-311`
  将 `robustness_threshold = 0.10` 直接当作 plotting band 的比例来源；
  这与 `REQ-29`、design 4.1 和 Python plot contract 冲突。
- `lwdid-r/tests/testthat/test-plot-diagnostics.R:348-358`
  要求 `GeomVline`，但 `REQ-31` 与 Python 只要求红色空心圆。

## 裁决

- `TC-10.2.23`: `confirmed test-oracle drift`
- `TC-10.2.27`: `confirmed test-oracle drift`
- `TC-10.2.2`: `confirmed implementation regression`

## 建议

1. `story-worker` 修 plot 时，应把折线层、`show_ci` 透传和 unconverged fallback
   视为实现问题继续修复。
2. `controller` / `correct-course` 在使用 `plot-diagnostics` 作为 gate 前，应先将
   `TC-10.2.23` 与 `TC-10.2.27` 降级为 oracle drift，或明确升格为
   R-only UX enhancement，而不能把它们包装成论文/Python parity 真值。
