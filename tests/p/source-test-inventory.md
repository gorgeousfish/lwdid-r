# Source Test Inventory

## 当前已识别的 R 侧测试资产

- `lwdid-r/tests/testthat/test-aggregate.R`
- `lwdid-r/tests/testthat/test-control-groups.R`
- `lwdid-r/tests/testthat/test-clustering-diagnostics.R`
- `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`
- `lwdid-r/tests/testthat/test-ipw.R`
- `lwdid-r/tests/testthat/test-ipwra.R`
- `lwdid-r/tests/testthat/test-psm.R`
- `lwdid-r/tests/testthat/test-randomization.R`
- `lwdid-r/tests/testthat/test-sensitivity-pre-period.R`
- `lwdid-r/tests/testthat/test-sensitivity-anticipation.R`
- `lwdid-r/tests/testthat/test-sensitivity-comprehensive.R`
- `lwdid-r/tests/testthat/test-sensitivity-comprehensive-common-controls-parity.R`
- `lwdid-r/tests/testthat/test-sensitivity-comprehensive-realdata-parity.R`
- `lwdid-r/tests/testthat/test-sensitivity-comprehensive-monte-carlo-parity.R`
- `lwdid-r/tests/testthat/test-sensitivity-comprehensive-hardening-parity.R`
- `lwdid-r/tests/testthat/test-sensitivity-comprehensive-fatal-inheritance-parity.R`
- `lwdid-r/tests/testthat/test-sensitivity-comprehensive-stata-precision-parity.R`
- `lwdid-r/tests/testthat/test-sensitivity-comprehensive-staggered-realdata-parity.R`
- `lwdid-r/tests/testthat/test-staggered-*.R`
- `lwdid-r/tests/testthat/test-stata-consistency.R`
- `lwdid-r/tests/testthat/test-lwdid-integration*.R`

## 当前已识别的 parity / comparator 资产

- `_automation/test-artifacts/parity/e8_03_common_timing_comparator.py`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-common-timing-comparator.json`
- `_automation/test-artifacts/parity/e8_03_smoking_realdata_comparator.py`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-comparator.json`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-comparator.md`
- `_automation/test-artifacts/parity/e8_03_smoking_stata_precision_comparator.py`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-precision.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-precision.md`
- `_automation/test-artifacts/parity/e8_03_smoking_stata_sensitivity_branches_comparator.py`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.md`
- `_automation/test-artifacts/parity/20260324-theory-parity-e8-03-pre-period-npre1-manual-oracle.json`
- `_automation/test-artifacts/parity/20260324-theory-parity-e8-03-pre-period-npre1-manual-oracle.md`
- `_automation/test-artifacts/parity/e8_03_comprehensive_numeric_audit.py`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-comprehensive-numeric-audit.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-comprehensive-numeric-audit.md`
- `_automation/test-artifacts/parity/e8_03_smoking_controls_comparator.py`
- `_automation/test-artifacts/parity/e8_03_smoking_frozen_controls_fixture.csv`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-frozen-controls-fixture.md`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-controls-comparator.json`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-smoking-controls-comparator.md`
- `_automation/test-artifacts/parity/e8_03_castle_realdata_comparator.py`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-castle-comparator.json`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-castle-comparator.md`
- `_automation/test-artifacts/parity/e8_03_monte_carlo_common_timing_comparator.py`
- `_automation/test-artifacts/parity/e8_03_monte_carlo_common_timing_fixture.csv`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-monte-carlo-common-timing.json`
- `_automation/test-artifacts/parity/20260323-qa-parity-e8-03-monte-carlo-common-timing.md`
- `_automation/test-artifacts/parity/e8_03_fatal_inheritance_audit.R`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-fatal-inheritance-regression.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-03-fatal-inheritance-regression.md`
- `_automation/test-artifacts/parity/e8_04_clustering_task1_oracle.py`
- `_automation/test-artifacts/parity/e8_04_clustering_hierarchical_fixture.csv`
- `_automation/test-artifacts/parity/e8_04_clustering_small_cluster_fixture.csv`
- `_automation/test-artifacts/parity/e8_04_clustering_within_cluster_variation_fixture.csv`
- `_automation/test-artifacts/parity/e8_04_clustering_never_treated_fixture.csv`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-clustering-task1-python-oracle.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-clustering-task1-python-oracle.md`
- `_automation/test-artifacts/parity/e8_04_clustering_task23_contract_oracle.py`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-clustering-task23-contract.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-clustering-task23-contract.md`
- `_automation/test-artifacts/parity/e8_04_clustering_layer3_state_policy_oracle.py`
- `_automation/test-artifacts/parity/e8_04_clustering_layer3_state_policy_fixture.csv`
- `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer3-state-policy.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer3-state-policy.md`
- `_automation/test-artifacts/parity/e8_04_clustering_layer3_few_cluster_oracle.py`
- `_automation/test-artifacts/parity/e8_04_clustering_layer3_few_cluster_fixture.csv`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer3-few-cluster.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer3-few-cluster.md`
- `_automation/test-artifacts/parity/e8_04_clustering_layer3_hierarchical_oracle.py`
- `_automation/test-artifacts/parity/e8_04_clustering_layer3_hierarchical_fixture.csv`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer3-hierarchical.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer3-hierarchical.md`
- `_automation/test-artifacts/parity/e8_04_clustering_layer3_smoking_common_timing_oracle.py`
- `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer3-smoking-common-timing.json`
- `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer3-smoking-common-timing.md`
- `_automation/test-artifacts/parity/e8_04_clustering_layer4_monte_carlo_contract.yaml`
- `_automation/test-artifacts/parity/e8_04_package_hardening_regression.R`
- `_automation/test-artifacts/parity/e8_04_layer5_release_regression.R`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer5-release-regression.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer5-release-regression.md`
- `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer4-coverage-large-g.json`
- `_automation/test-artifacts/parity/20260324-story-worker-e8-04-layer4-coverage-large-g.md`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer4-standard-scenarios.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer4-standard-scenarios.md`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer4-wild-bootstrap-scenarios.json`
- `_automation/test-artifacts/parity/20260324-qa-parity-e8-04-layer4-wild-bootstrap-scenarios.md`

## E8-04 clustering 当前来源映射

- Python archived source:
  `lwdid-py_v0.2.3/tests/_archived_dev_tests/root_level/test_clustering_diagnostics.py`
- Python implementation source:
  `lwdid-py_v0.2.3/src/lwdid/clustering_diagnostics.py`
- Shared fixtures:
  `e8_04_clustering_hierarchical_fixture.csv`,
  `e8_04_clustering_within_cluster_variation_fixture.csv`,
  `e8_04_clustering_never_treated_fixture.csv`,
  `e8_04_clustering_small_cluster_fixture.csv`,
  `e8_04_clustering_layer3_state_policy_fixture.csv`,
  `e8_04_clustering_layer3_few_cluster_fixture.csv`
- R regressions:
  `lwdid-r/tests/testthat/test-clustering-diagnostics.R`,
  `lwdid-r/tests/testthat/test-clustering-diagnostics-parity.R`
- 当前 coverage:
  `Task E8-04.1` 的 `.determine_cluster_level()` 与 `.analyze_cluster_var()`
  已具备 shared-fixture exact / numeric parity；`Task E8-04.2` 的
  `.detect_treatment_variation_level()` / internal recommender 也已有 helper-level
  regression；`Task E8-04.3` 则新增 public-contract oracle，冻结了
  small-cluster diagnosis / recommender anchors 与
  `check_clustering_consistency()` 的 local-scope + `dropna=True` 语义。
  `2026-03-24 08:51 CST qa-parity` 又通过
  `20260324-story-worker-e8-04-layer3-state-policy.json` +
  `test-clustering-diagnostics-parity.R` 锁定了首个 Layer 3 state-policy
  public-workflow parity，覆盖 `diagnose_clustering()` /
  `recommend_clustering()` / `check_clustering_consistency()` 的 end-to-end
  regression。`2026-03-24 10:07 CST qa-parity` 又通过
  `20260324-qa-parity-e8-04-layer3-few-cluster.json` +
  `test-clustering-diagnostics-parity.R` 锁定了 archived Python
  `TestFewClusters` 的 end-to-end public-workflow parity，并把
  `<20 clusters` 的 diagnose warning string 一并收口到 Python contract。
  `2026-03-24 12:33 CST qa-parity` 又通过
  `20260324-qa-parity-e8-04-layer3-hierarchical.json` +
  `test-clustering-diagnostics-parity.R` 锁定了 archived Python
  `test_clustering_empirical.py` 的 hierarchical panel public workflow：
  当前 `diagnose_clustering()` 会检测 `treatment_variation_level = "region"`，
  但 public `recommend_clustering()` 仍以 20 个 `state` clusters 的更强
  reliability / no-WCB contract 推荐 `state`，而
  `check_clustering_consistency(cluster_var = "state")` 保持 `0%` inconsistency。
  `2026-03-24 12:42 CST story-worker` 又通过
  `20260324-story-worker-e8-04-layer3-smoking-common-timing.json` +
  `test-clustering-diagnostics-parity.R` 锁定了 `smoking` common-timing
  built-in empirical public workflow：当前 `diagnose_clustering()`、
  `recommend_clustering()` 与 `check_clustering_consistency(cluster_var = "state")`
  都与 Python 对齐为 `treatment_variation_level = "state"`、
  `recommended_var = "state"`、`n_clusters = 39`、`0% inconsistency`、
  `use_wild_bootstrap = FALSE`。
- 当前剩余 gap:
  当前 synthetic Layer 2、Layer 3 state-policy、Layer 3 few-cluster、
  Layer 3 hierarchical 与 Layer 3 smoking common-timing regression 已不再缺
  fixture / oracle；Layer 4 的 standard cluster-robust path 也已由
  `20260324-story-worker-e8-04-layer4-coverage-large-g.json/.md`、
  `20260324-qa-parity-e8-04-layer4-standard-scenarios.json/.md` 与
  `20260324-qa-parity-e8-04-layer4-wild-bootstrap-scenarios.json/.md`
  固定为 executable regression。`story-E8-04` 当前已无剩余 parity gap；但
  fresh `devtools::check(--no-manual)` 说明 story 仍被 package-hardening blocker
  卡住。`20260324-qa-parity-e8-04-layer5-release-regression.json/.md` 已确认
  `TC-10.6.18` / `TC-10.6.34` 不再是 live blocker，而
  `20260324-qa-parity-e8-04-layer5-tests-rebucket.json/.md` 又进一步确认：
  `transform-detrend` / `utils` / `validate` 三条 stale tests contract 当前也已转绿；
  当前 canonical full-check tests bucket 已改写为
  `test-staggered-integration.R:383/531` 与
  `test-staggered-numerical.R:181` 三条 check-time-only aggregation contract，
  同时 remaining non-tests bucket 继续停留在 `R code for possible problems`
  与 `Rd files / Rd cross-references / missing documentation / codoc / usage`
  hygiene。下一步不应直接切换到 `story-E8-05`，而应先处理这些 closure-ready
  blocker。

## 当前需要系统补齐的来源映射

1. Python tests -> R tests 的一一对标表
2. 论文数值示例 -> comparator / oracle 列表
3. `story-E8-03` 当前已有 synthetic、`smoking` real-data、`castle`
   real-data 与 common-timing Monte Carlo 四条共享证据路径；Python 侧仍缺直接的 comprehensive
   test fixture 映射表
4. 发布前 regression / Stata / FATAL 收口表
   - `castle` unit-level `X_i` controls bug-aware comparator 已冻结
   - `smoking` common-timing / no-controls transformation path 已新增 Stata-backed precision oracle
   - `smoking` comprehensive 的 `pre_period` / `no_anticipation` Stata oracle 已新增：`pre_period n_pre=2:19` 与 `exclude=0:3` 当前通过；`n_pre=1` 虽在 Stata 端返回 `r(2001)`，但 `20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json/.md` 已把 theory-parity 的 paper-backed manual oracle 接入同一 QA contract，并由 `test-sensitivity-comprehensive-stata-precision-parity.R` 对 `ATT/SE/p/CI` 做 executable regression，因此 `8.4` 已收口
   - `8.6` 已新增 comprehensive hardening audit：覆盖 `smoking`、`castle` 与 `smoking-frozen-controls` 三条 `type="all"` real-data 路径
   - `smoking` common-timing shared frozen `X_i` controls comparator 已冻结
   - `8.5` 已新增 fatal inheritance audit：`FATAL-001/002/004` 当前具备 executable regression
   - 后续主要缺口已不再是 controls fixture、Task `8.4`、generic FATAL inheritance 或 `story-E8-04` 的 Layer 4 Monte Carlo；`story-E8-03` 与 `story-E8-04` 的 parity 证据当前都已具备 executable coverage。系统层面的剩余缺口转为 Python tests -> R tests 的来源映射表，以及 `E8-05` 的 helper-public API parity 入口

## 使用要求

- qa-parity 每轮若新增 comparator、oracle、fixture 或 tolerance 证据，应同步更新本文件。
