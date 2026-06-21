# story-E8-04 边界状态收敛

- 时间：`2026-03-24 12:19 CST`
- 角色：`correct-course`
- phase：`implementation`
- active_story：`story-E8-04`

## 本轮目标

确认 `story-E8-04` 当前是否仍真实阻塞在 `Task 4.4/4.8/5.2`，还是这些边界矩阵已经由
现有代码与测试收口，只剩 Layer 3 / Layer 4 parity gap。

## 新鲜证据

1. fresh rerun
   `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="clustering-diagnostics|clustering-diagnostics-parity", reporter="summary")'`
   当前通过。
2. `lwdid-r/tests/testthat/test-clustering-diagnostics.R` 已直接覆盖：
   - `E8-04.4.4`：`4% / 5% / 0% / 100%` inconsistency 边界；
   - `E8-04.4.8`：nested / `G=1` / `level="same"` 的 `is_valid` 三条件矩阵；
   - `E8-04.5.2`：单候选、空候选报错、缺失变量报错、`CV=0`、all-never-treated。
3. `lwdid-r/R/clustering_diagnostics.R` 当前合同保持为：
   - `pct_inconsistent < 5`
   - `is_valid = !is_nested_in_unit && n_clusters >= 2L && level != "lower"`
   - small-cluster warning string 已与 Python-backed few-cluster oracle 对齐。

## 结论

- `confirmed`：当前不存在新的 clustering 业务 drift。
- `confirmed`：`Task 4.4`、`Task 4.8`、`Task 5.2` 已具备 executable evidence，
  不应继续写入 blocker。
- `remaining-work`：`story-E8-04` 当前只剩 `Task 5.3` 的 hierarchical /
  built-in empirical / common-timing supplemental Layer 3，以及 Layer 4
  Monte Carlo acceptance-band parity。

## 建议

1. `qa-parity` 下一轮直接抽取 hierarchical / built-in empirical /
   common-timing public-workflow oracle，不再重复 few-cluster。
2. `controller` / `correct-course` 后续状态收敛不得再把 `4.4/4.8/5.2`
   写回 blocker。
