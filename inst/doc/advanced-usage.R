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
#   controls = c("police", "income", "poverty"),
#   aggregate = "cohort"
# )

## ----ipwra, eval=FALSE--------------------------------------------------------
# result_ipwra <- lwdid(
#   data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#   gvar = "gvar", rolling = "demean", estimator = "ipwra",
#   controls = c("police", "income", "poverty"),
#   aggregate = "cohort"
# )

## ----psm, eval=FALSE----------------------------------------------------------
# result_psm <- lwdid(
#   data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#   gvar = "gvar", rolling = "demean", estimator = "psm",
#   controls = c("police", "income", "poverty"),
#   aggregate = "cohort"
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
#   wcb_reps = 999
# )
# summary(result_wcb)

## ----ri, eval=FALSE-----------------------------------------------------------
# result_ri <- lwdid(
#   data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
#   d = "d", post = "post", ri = TRUE, rireps = 500
# )
# cat("RI p-value:", result_ri$ri_pvalue, "\n")

## ----parallel, eval=FALSE-----------------------------------------------------
# result_par <- lwdid(
#   data = castle, y = "lhomicide", ivar = "sid", tvar = "year",
#   gvar = "gvar", rolling = "demean", aggregate = "cohort",
#   parallel = TRUE, n_cores = 4
# )

## ----parallel-setup, eval=FALSE-----------------------------------------------
# setup_parallel(n_workers = 4, strategy = "multisession")
# # All subsequent lwdid() calls use parallel backend
# result <- lwdid(...)
# # Reset when done
# setup_parallel(strategy = "sequential")

## ----export, eval=FALSE-------------------------------------------------------
# to_latex(res_dm, file = "results_table.tex")
# to_csv(res_dm, file = "results.csv")

