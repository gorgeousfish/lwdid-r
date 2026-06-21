## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----load---------------------------------------------------------------------
library(lwdid)
data(smoking)
data(castle)

## ----ipw, eval=FALSE----------------------------------------------------------
# result_ipw <- lwdid(
#   data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#   gvar = "gvar", rolling = "demean", estimator = "ipw",
#   ps_controls = c("police", "income", "poverty"),
#   aggregate = "cohort", control_group = "never_treated"
# )

## ----ipwra, eval=FALSE--------------------------------------------------------
# result_ipwra <- lwdid(
#   data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#   gvar = "gvar", rolling = "demean", estimator = "ipwra",
#   controls = c("police", "income", "poverty"),
#   ps_controls = c("police", "income", "poverty"),
#   aggregate = "cohort", control_group = "never_treated"
# )

## ----psm, eval=FALSE----------------------------------------------------------
# result_psm <- lwdid(
#   data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#   gvar = "gvar", rolling = "demean", estimator = "psm",
#   ps_controls = c("police", "income", "poverty"),
#   aggregate = "cohort", control_group = "never_treated"
# )

## ----transforms---------------------------------------------------------------
res_dm <- lwdid(
  data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
  d = "d", post = "post", rolling = "demean"
)

res_dt <- lwdid(
  data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
  d = "d", post = "post", rolling = "detrend"
)

cat("Demean ATT:", round(res_dm$att, 4), "\n")
cat("Detrend ATT:", round(res_dt$att, 4), "\n")

## ----vce-example--------------------------------------------------------------
# HC1 for smoking data (single treated unit makes HC3 unstable)
result_hc1 <- suppressWarnings(lwdid(
  data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
  d = "d", post = "post", vce = "hc1"
))
summary(result_hc1)

## ----compare-basic, eval=FALSE------------------------------------------------
# # Estimate with different estimators
# res_ra <- lwdid(data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#                 gvar = "gvar", rolling = "demean", estimator = "ra",
#                 aggregate = "overall")
# res_ipw <- lwdid(data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#                  gvar = "gvar", rolling = "demean", estimator = "ipw",
#                  ps_controls = c("police", "income", "poverty"),
#                  aggregate = "overall")
# res_ipwra <- lwdid(data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#                    gvar = "gvar", rolling = "demean", estimator = "ipwra",
#                    controls = c("police", "income", "poverty"),
#                    ps_controls = c("police", "income", "poverty"),
#                    aggregate = "overall")
# 
# # Compare specifications side by side
# # Output: table with ATT, SE, CI, significance stars, model metadata
# compare(RA = res_ra, IPW = res_ipw, IPWRA = res_ipwra)

## ----compare-transform, eval=FALSE--------------------------------------------
# # Demean vs Detrend comparison
# res_demean <- lwdid(data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#                     gvar = "gvar", rolling = "demean", estimator = "ra",
#                     aggregate = "overall")
# res_detrend <- lwdid(data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#                      gvar = "gvar", rolling = "detrend", estimator = "ra",
#                      aggregate = "overall")
# 
# # Similar estimates suggest robustness to trend specification
# compare(Demean = res_demean, Detrend = res_detrend)

## ----compare-effects, eval=FALSE----------------------------------------------
# # Effects-level comparison shows cohort-by-cohort differences
# compare(RA = res_ra, IPWRA = res_ipwra, type = "effects")
# 
# # Output columns include:
# # - Model name, ATT estimate, standard error
# # - 95% confidence interval [ci_lower, ci_upper]
# # - Significance stars (*** p<0.01, ** p<0.05, * p<0.1)
# # - Model metadata: estimator, transformation, N observations

## ----sensitivity--------------------------------------------------------------
sens <- suppressWarnings(lwdid_sensitivity(
  data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
  d = "d", post = "post", type = "all", verbose = FALSE
))
print(sens)

## ----pre-period-sens, eval=FALSE----------------------------------------------
# pre_sens <- lwdid_sensitivity_pre_period(
#   data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
#   d = "d", post = "post"
# )

## ----pt-test------------------------------------------------------------------
pt <- lwdid_test_parallel_trends(
  data = castle, y = "lhomicide", ivar = "sid",
  tvar = "year", gvar = "gvar"
)
print(pt)

## ----recommend----------------------------------------------------------------
rec <- lwdid_recommend_transformation(
  data = castle, y = "lhomicide", ivar = "sid",
  tvar = "year", gvar = "gvar"
)
print(rec)

## ----diagnose-----------------------------------------------------------------
diag <- suppressWarnings(lwdid_diagnose(
  data = castle, y = "lhomicide", ivar = "sid",
  tvar = "year", gvar = "gvar"
))
summary(diag)

## ----aggregate, eval=FALSE----------------------------------------------------
# agg_result <- aggregate_to_panel(
#   data = survey_data,
#   unit_var = "state",
#   time_var = "year",
#   outcome_var = "income",
#   weight_var = "survey_weight",
#   treatment_var = "treated",
#   min_cell_size = 30
# )
# panel <- as.data.frame(agg_result)

## ----wcb, eval=FALSE----------------------------------------------------------
# result_wcb <- lwdid(
#   data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
#   d = "d", post = "post", rolling = "demean",
#   vce = "bootstrap", cluster_var = "state",
#   wcb_reps = 999, wcb_seed = 42,
#   use_fwildclusterboot = FALSE
# )
# summary(result_wcb)

## ----ri, eval=FALSE-----------------------------------------------------------
# result_ri <- lwdid(
#   data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
#   d = "d", post = "post", ri = TRUE, ri_method = "permutation",
#   rireps = 500, seed = 42
# )
# cat("RI p-value:", result_ri$ri_pvalue, "\n")

## ----parallel, eval=FALSE-----------------------------------------------------
# result_par <- lwdid(
#   data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#   gvar = "gvar", rolling = "demean", aggregate = "cohort",
#   control_group = "never_treated", parallel = TRUE, n_cores = 4
# )

## ----parallel-setup, eval=FALSE-----------------------------------------------
# setup_parallel(n_workers = 4, strategy = "multisession")
# # All subsequent lwdid() calls use parallel backend
# result <- lwdid(...)
# # Reset when done
# setup_parallel(strategy = "sequential")

## ----tidy-example, eval=FALSE-------------------------------------------------
# library(lwdid)
# result <- lwdid(data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#                 gvar = "gvar", rolling = "demean", estimator = "ra",
#                 aggregate = "overall")
# 
# # Overall ATT (default)
# # Returns: data.frame with columns term, estimate, std.error, statistic,
# #          p.value, conf.low, conf.high
# tidy(result)
# 
# # Cohort-level effects
# # Returns one row per cohort with cohort-specific ATT
# tidy(result, type = "effects")
# 
# # Event-time aggregation
# # Returns one row per relative event-time period
# tidy(result, type = "event_time")

## ----glance-example, eval=FALSE-----------------------------------------------
# # Returns: nobs, n_treated, n_control, n_cohorts, n_periods,
# #          estimator, transformation, vce_type, r.squared (if applicable)
# glance(result)

## ----fixest-style, eval=FALSE-------------------------------------------------
# # Standard errors
# se(result)
# 
# # t-statistics
# tstat(result)
# 
# # p-values
# pvalue(result)
# 
# # Full coefficient table (estimate, se, t-stat, p-value)
# coeftable(result)

## ----modelsummary-example, eval=FALSE-----------------------------------------
# library(modelsummary)
# 
# # Compare multiple specifications in one table
# res_ra <- lwdid(data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#                 gvar = "gvar", rolling = "demean", estimator = "ra",
#                 aggregate = "overall")
# res_ipwra <- lwdid(data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#                    gvar = "gvar", rolling = "demean", estimator = "ipwra",
#                    controls = c("police", "income", "poverty"),
#                    ps_controls = c("police", "income", "poverty"),
#                    aggregate = "overall")
# 
# # Produces a formatted comparison table with estimates, SEs, and fit stats
# modelsummary(list(RA = res_ra, IPWRA = res_ipwra))
# 
# # Export to LaTeX or Word
# modelsummary(list(RA = res_ra, IPWRA = res_ipwra),
#              output = "results.tex")
# modelsummary(list(RA = res_ra, IPWRA = res_ipwra),
#              output = "results.docx")

## ----export, eval=FALSE-------------------------------------------------------
# to_latex(res_dm, file = "results_table.tex")
# to_csv(res_dm, file = "results.csv")

