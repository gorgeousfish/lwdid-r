# E8-04 Layer 4 wild-bootstrap source contract

## 本轮结论

- `confirmed`: `wild_bootstrap_coverage_small_g` 与
  `wild_bootstrap_relative_gap` 这两个剩余 blocker 不是简单的
  “`n_bootstrap = 199` 随机抽样”合同。Python 在 `G = 10` 且默认
  `weight_type = "rademacher"` 时，会自动启用 full enumeration，实际枚举
  `2^10 = 1024` 个权重组合。
- `confirmed`: archived Python tests 调用的 wild bootstrap 默认保持
  `impose_null = TRUE`，并用 `percentile_t` 方式构造置信区间；p-value 采用
  `P(|t*| >= |t_orig|)` 的包含式规则，而不是严格 `>`。
- `confirmed`: 这些都是 Python-backed QA contract 的实现细节，不是
  `Docs/lw2026.md` 的论文定理阈值；因此不更新 `Docs/Python包bug列表.md`。

## 直接 source anchors

### Archived Python tests 如何调用 WCB

- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:165-170`
  在 `wild_bootstrap_coverage_small_g` 场景调用：
  `wild_cluster_bootstrap(data, y_transformed='Y', d='D', cluster_var='cluster', n_bootstrap=199, alpha=0.05)`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:215-219`
  在 `wild_bootstrap_relative_gap` 场景调用：
  `wild_cluster_bootstrap(data, ..., n_bootstrap=199)`。

### Python 实际执行的 WCB 机制

- `lwdid-py_v0.2.3/src/lwdid/inference/wild_bootstrap.py:726-788`
  定义默认参数：`weight_type='rademacher'`、`alpha=0.05`、
  `impose_null=True`、`full_enumeration=None`，并明确写明
  `G <= 12` 且 `weight_type='rademacher'` 时会自动启用 full enumeration。
- `lwdid-py_v0.2.3/src/lwdid/inference/wild_bootstrap.py:860-940`
  实际把 `full_enumeration` 解析为
  `G <= 12 and weight_type == 'rademacher'`；在 `G = 10` 的两个剩余 scenario 中，
  `actual_n_bootstrap = len(all_weights) = 2^10 = 1024`，因此请求的
  `n_bootstrap = 199` 不再是最终执行次数。
- `lwdid-py_v0.2.3/src/lwdid/inference/wild_bootstrap.py:1019-1054`
  定义 p-value 与 CI：
  - p-value: `mean(abs(t_stats_valid) >= abs(t_stat_original))`
  - `impose_null = TRUE` 时 CI 使用对称 `percentile_t`
  - 返回结果中的 `n_bootstrap` 是 `actual_n_bootstrap`

### 当前 R 侧状态

- `rg "wild_cluster|wild bootstrap|full_enumeration|impose_null|percentile_t" lwdid-r/R`
  当前除 clustering recommendation 文案外，没有对应的 WCB 实现入口；
  因此 `story-E8-04` 的剩余 gap 是明确的实现缺口，而不是来源不清。

## 对下一轮的裁决

1. 后续 R parity 若要实现这两个 wild-bootstrap scenario，必须优先复现 Python 的
   `effective contract`：`weight_type = rademacher`、
   `impose_null = TRUE`、`percentile_t`、`|t*| >= |t_orig|`，以及
   `G = 10` 时的 `2^10 = 1024` full enumeration。
2. 不得把剩余 blocker 简化成“199 次随机 bootstrap”；那会引入额外 Monte Carlo
   噪声，并偏离 archived Python tests 的真实执行路径。
3. 由于这些细节来自 Python WCB implementation，而不是论文直接定理，发布层文案
   只能写成 “Python-backed QA contract passed/failed”，不能写成
   “paper theorem satisfied/violated”。
