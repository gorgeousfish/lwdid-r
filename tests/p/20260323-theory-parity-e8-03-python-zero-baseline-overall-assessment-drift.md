# story-E8-03 Python comprehensive sensitivity 零基线总体判断漂移证据

## 结论

- `confirmed`: `lwdid-py_v0.2.3/src/lwdid/sensitivity.py` 在 transformation
  comparison 中对 `recommendations` 使用了 `abs(demean_att) > 1e-10`
  的零基线保护，但在 `overall_assessment` 的 `issues` 判定中绕开了该保护。
- `confirmed`: 当 `demean_att` 接近零而 `difference > 0.25 * |demean_att|` 时，
  Python 会同时返回
  `Caution: transformation sensitivity detected. See recommendations.`
  和 `No major robustness concerns identified.`，形成同一结果对象内部的
  自相矛盾。
- `not-a-paper-rule`: 论文支持比较 demean 与 detrend 的 transformed ATT，
  但并未给出“25% 相对差异”这一 heuristic，更没有把“25% of zero”当作可解释
  的理论对象。因此近零基线时不应把 Python 的该字符串输出当成真值。

## 真值锚点

### 1. 论文中的 invariant object

- `Docs/lw2026.md:198-236` 只把 demean 与 detrend 比较建立在 transformed
  outcomes 与 ATT 估计之上：
  - demean: Procedure 2.1 的 transformed outcome；
  - detrend: Procedure 3.1 的 $\ddot{Y}_{it}$ 与 $\hat{\tau}_{DT}$。
- `Docs/lw2025.md:463-517` 说明 detrending 的作用是当存在异质线性趋势时，
  先消除 trend component，再在 transformed outcome 上做 ATT 识别。

结论：论文支持比较两种 transformed ATT 是否差异明显，但没有定义一个
 “相对零基线”的理论比率。因此近零基线只能视为 heuristic 的数值边界问题，
 不能当作论文结论。

### 2. Python 自身的一致性基准

- `lwdid-py_v0.2.3/src/lwdid/sensitivity.py:1346-1351` 中
  `_compute_sensitivity_ratio()` 已明确把 near-zero baseline 当成单独的数值边界：
  - 若 `abs(baseline_att) > 1e-10`，返回相对比率；
  - 否则仅在 `att_range > 1e-10` 时返回 `Inf`，否则返回 `0.0`。

这说明 Python 包内部本来就承认“接近零的 baseline 不能直接做普通相对比例比较”。

## 源码证据

### Python 当前行为

- `lwdid-py_v0.2.3/src/lwdid/sensitivity.py:2279-2294`
  在 transformation comparison 中构造 `difference`，并只在
  `abs(result_demean.att) > 1e-10` 时计算 `rel_diff` 和添加 recommendation。
- `lwdid-py_v0.2.3/src/lwdid/sensitivity.py:2340-2347`
  在 `overall_assessment` 中却直接使用
  `difference > 0.25 * abs(demean_att)` 来加入
  `transformation sensitivity` issue，没有复用上面的 guard。

### R 当前行为

- `lwdid-r/R/sensitivity.R:1453-1465`
  先生成 `difference` 与 guarded `rel_diff`。
- `lwdid-r/R/sensitivity.R:1572-1583`
  只有在 `rel_diff` 非 `NA` 且 `rel_diff > 0.25` 时，才把它记为
  `变换方法敏感性`。

结论：R 当前实现没有跟随 Python 的该漂移。

## 最小复现

以下脚本在不依赖真实估计器数值的前提下，直接复现 Python 的综合结果内部矛盾：

```python
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path("lwdid-py_v0.2.3/src").resolve()))
from lwdid.sensitivity import sensitivity_analysis

df = pd.DataFrame({
    "id": [1, 1, 2, 2],
    "time": [1, 2, 1, 2],
    "y": [0.0, 0.0, 0.0, 0.0],
    "d": [0, 1, 0, 0],
    "post": [0, 1, 0, 1],
})

def fake_lwdid(**kwargs):
    if kwargs["rolling"] == "demean":
        return SimpleNamespace(att=1e-12, se_att=0.1)
    if kwargs["rolling"] == "detrend":
        return SimpleNamespace(att=2e-12, se_att=0.1)
    raise AssertionError(kwargs["rolling"])

with patch("lwdid.sensitivity._validate_robustness_inputs", return_value=None), \
     patch("lwdid.sensitivity._get_max_pre_periods", return_value=2), \
     patch("lwdid.core.lwdid", side_effect=fake_lwdid):
    result = sensitivity_analysis(
        data=df,
        y="y",
        ivar="id",
        tvar="time",
        d="d",
        post="post",
        analyses=["transformation"],
        verbose=False,
    )

print(result.overall_assessment)
print(result.recommendations)
```

复现结果：

- `overall_assessment`:
  `Caution: transformation sensitivity detected. See recommendations.`
- `recommendations`:
  `['No major robustness concerns identified.']`

## 裁决

- `Python bug status`: `confirmed`
- `影响范围`:
  - 近零 `demean_att` 的 comprehensive sensitivity summary；
  - 后续 parity 若用 exact string 对比 `overall_assessment`，会把 Python 的
    自相矛盾输出误判为真值。
- `R 是否跟随`: 否。R 应继续以 guarded `rel_diff` 为准。
- `parity 处理建议`:
  - 近零基线场景下，不要求 `overall_assessment` 与 Python exact string parity；
  - 应比较 transformed ATT、difference、以及 guard 后的 issue 判定是否自洽；
  - 如需与 Python 对齐，只能以 bug ledger 的 waived 方式说明，不能倒逼 R 回退。
