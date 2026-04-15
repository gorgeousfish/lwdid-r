# E8-04 AC-17 fresh check 审计（story-worker）

- 时间：2026-03-24 17:05 CST
- 角色：story-worker
- active_story：story-E8-04
- phase：implementation

## 本轮结论

`story-E8-04` 的 parity / clustering 回归仍保持通过，但 `AC-17` 还没有收口。
本轮 fresh `devtools::check(pkg=".../lwdid-r", args="--no-manual", error_on="never")`
显示：旧的 `compute_inference.Rd` parse、`aggregate_to_overall.Rd` parse、
`test-sensitivity-pre-period.R` parse、`namespace nobs` load failure 当前都未再出现；
真实 blocker 已切换为更上层的 package hardening 清单。

## 本轮新增工件

- `_automation/test-artifacts/parity/e8_04_compute_inference_rd_regression.R`
  新增 `compute_inference.Rd` regression guard：
  1. 禁止再次出现会触发 Rd 注释吞噬的 `\\%` 双转义序列
  2. 要求 `tools::Rd2txt()` 可无 warning 转换

## fresh check 剩余 blocker

### WARNING

1. `code files for non-ASCII characters`
   - `R/clustering_diagnostics.R`
   - `R/plot_diagnostics.R`
   - `R/plot_event_study.R`
   - `R/sensitivity.R`
   - `R/staggered.R`
2. `dependencies in R code`
   - 未声明 `MASS`、`gridExtra`、`patchwork`、`scales`
3. `Rd files`
   - `dot-escape_latex.Rd:18`
   - `dot-filter_excluding_periods.Rd:26`
   - `grapes-or-or-grapes.Rd:3`
   - `precompute_transforms.Rd:57`
4. `Rd cross-references`
   - `aggregate_to_event_time.Rd`
   - `apply_precomputed_transform.Rd`
   - `dot-filter_excluding_periods.Rd`
   - `estimate_outcome_model.Rd`
   - `estimate_propensity_score.Rd`
   - `get_unit_level_gvar.Rd`
5. `missing documentation entries`
   - `to_latex_comparison`
6. `code/documentation mismatches`
   - `dot-new_lwdid_sensitivity_comprehensive.Rd`
7. `Rd usage sections`
   - `dot-estimate_common_timing.Rd`
   - `lwdid-conditions.Rd`
   - `lwdid.Rd`
   - `new_lwdid_result.Rd`
8. `unstated dependencies in tests`
   - `devtools`
   - `jsonlite`
   - `yaml`

### NOTE

1. `Packages suggested but not available for checking`
   - `fwildclusterboot`
   - `openxlsx`
2. `Vignette dependency required without any vignettes`
   - `knitr`
3. `unable to verify current time`
4. `R code for possible problems`
   - 以未导入 `stats` / `utils` 函数与 data.table NSE 变量为主

## 对下一轮的约束

- 不要再把 `compute_inference.Rd` parse / `test-sensitivity-pre-period.R` parse /
  `namespace nobs` 载入失败写回 canonical blocker。
- 下一轮更合适的 bounded package 是：
  1. 先处理 roxygen unresolved-link / Rd usage / codoc 集群
  2. 再收敛 imports 与 non-ASCII
