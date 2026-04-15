# story-E8-04 release / state convergence

## 时间

- `2026-03-24 18:18:53 CST`

## fresh rerun

- `Rscript /Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e8_04_layer5_release_regression.R`
  当前通过。
- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="visualization-export", reporter="summary")'`
  当前通过；仅保留 `TC-10.6.51/52/53` 的 warning，不再有 failure。
- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="vce-integration", reporter="summary")'`
  当前通过。
- `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="estimate-ra-common", reporter="summary")'`
  初始仅剩 `T-32` 失败，根因是测试仍要求 `estimate_ra_common()` 返回 19 个字段；而当前 canonical contract 已由
  `test-vce-return-values.R` 固定为 21 个字段（新增 `df_resid` /
  `df_inference`）。

## 最小修复

- 已将 `lwdid-r/tests/testthat/test-estimate-ra-common.R` 的 `T-32` 期望从 19
  字段同步到 21 字段，避免 test contract 长期落后于代码事实。
- 修正后 fresh rerun
  `Rscript -e 'devtools::test(pkg="/Users/cxy/Desktop/lwdid_r/lwdid-r", filter="estimate-ra-common", reporter="summary")'`
  当前通过。

## 裁决

- `TC-10.6.18` 与 `TC-10.6.34` 不再是 live blocker；旧的
  `visualization-export` release-contract 叙述已过期。
- 当前 `story-E8-04` 仍不能写成 closure-ready，但 blocker 需要回到 fresh full
  `devtools::check(--no-manual)` 重新确认；截至 `2026-03-24 17:28 CST`
  的最近一次 full check，remaining bucket 主要是 `R code globals note` 与
  `Rd/codoc/missing documentation/usage` hygiene。
- 因此本轮纠偏同时完成：
  1. 一个真实 test contract drift 修正；
  2. 一个 release-state drift 证据收敛。
