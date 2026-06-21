# E8-04 Layer 4 Monte Carlo contract

## 本轮结论

- `confirmed`: `story-E8-04` 的唯一剩余 blocker 已可压缩成
  `test_clustering_monte_carlo.py` 的 9 个 band-based scenarios；下一轮不需要再
  重新解释 Layer 4 来源。
- `confirmed`: `_automation/test-artifacts/parity/e8_04_clustering_layer4_monte_carlo_contract.yaml`
  已把各 scenario 的 DGP、`seed`、`n_simulations`、`G`、
  `obs_per_cluster`、`n_bootstrap` 与 acceptance bands 固化成 machine-readable
  contract。
- `confirmed`: `Docs/lw2026.md:641-643` 只提供“处理在更高层级分配时可在更高层级
  聚类”的论文真值；Layer 4 的 coverage / size / power 数值阈值并非 paper-only
  truth，而是 Python archived QA contract。
- `confirmed`: 本轮未发现新的 Python-paper drift，因此不更新
  `Docs/Python包bug列表.md`。

## 直接 source anchors

### 论文边界

- `Docs/lw2026.md:641-643`
  只说明：在 `lwdid` 变换后的横截面回归中，若政策只在高于观测单位的层级变化，
  且 treated / control clusters 足够多，则可在该更高层级聚类。
- `Docs/lw2026.md:286-405`
  的 Monte Carlo 章节讨论的是 detrending estimator 的一般有限样本表现与
  nominal 95% coverage 附近的论文叙述，不提供 clustering diagnostics 模块的
  `G = 50/20/10` coverage bands、size bands 或 WCB bands。

### Python QA contract

- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:29-64`
  固定 `G = 50`、`obs_per_cluster = 30` 的 coverage band
  `0.90 <= coverage_rate <= 0.99`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:66-97`
  固定 `G = 20`、`obs_per_cluster = 50` 的 coverage band
  `0.88 <= coverage_rate <= 0.99`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:99-134`
  固定 `G = 10` small-cluster baseline：`coverage_rate > 0.80`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:141-177`
  固定 `G = 10` wild bootstrap coverage：`coverage_rate > 0.85`，
  且 `n_bootstrap = 199`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:179-230`
  固定 relative-gap contract：
  `coverage_wild_rate >= coverage_std_rate - 0.10`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:240-274`
  固定 size-under-null band：`0.02 <= rejection_rate <= 0.10`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:276-310`
  固定 power-under-alternative band：`power > 0.50`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:320-352`
  固定 unbalanced-treatment robustness：`coverage_rate > 0.85`。
- `lwdid-py_v0.2.3/tests/_archived_dev_tests/test_clustering_monte_carlo.py:354-402`
  固定 heterogeneous-cluster-sizes robustness：`coverage_rate > 0.85`。

## 对下一轮的裁决

1. Layer 4 parity 应严格读取
   `_automation/test-artifacts/parity/e8_04_clustering_layer4_monte_carlo_contract.yaml`
   的 scenario matrix，而不是重新自拟 bands。
2. Monte Carlo parity 只能写成区间断言，不能写成 exact-number oracle。
3. 若后续 R 侧 Monte Carlo 结果偏离这些 bands，先判断是 R drift、实现差异，
   还是 Monte Carlo noise；不要直接把 Python archived test 的单次输出当真值。
4. 由于这些阈值本身不是论文定理，发布文档只能说“Python-backed QA target
   passed/failed”，不能说“paper theorem satisfied/violated”。
