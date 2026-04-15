## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
# # From GitHub
# # install.packages("devtools")
# devtools::install_github("gorgeousfish/lwdid-r")

## ----common-timing-basic------------------------------------------------------
library(lwdid)
data(smoking)

result <- lwdid(
  data = smoking,
  y = "lcigsale",
  ivar = "state",
  tvar = "year",
  d = "d",
  post = "post",
  rolling = "demean",
  estimator = "ra"
)

print(result)

## ----common-timing-detrend----------------------------------------------------
result_dt <- lwdid(
  data = smoking,
  y = "lcigsale",
  ivar = "state",
  tvar = "year",
  d = "d",
  post = "post",
  rolling = "detrend",
  estimator = "ra"
)

summary(result_dt)

## ----staggered-basic----------------------------------------------------------
data(castle)

result_stag <- lwdid(
  data = castle,
  y = "lhomicide",
  ivar = "sid",
  tvar = "year",
  gvar = "gvar",
  rolling = "demean",
  estimator = "ra",
  aggregate = "cohort"
)

print(result_stag)

## ----staggered-event-study, eval=requireNamespace("ggplot2", quietly=TRUE)----
result_event <- lwdid(
  data = castle,
  y = "lhomicide",
  ivar = "sid",
  tvar = "year",
  gvar = "gvar",
  rolling = "demean",
  estimator = "ra",
  aggregate = "event_time"
)

plot(result_event)

## ----controls, eval=FALSE-----------------------------------------------------
# result_ctrl <- lwdid(
#   data = castle,
#   y = "lhomicide",
#   ivar = "sid",
#   tvar = "year",
#   gvar = "gvar",
#   rolling = "demean",
#   estimator = "ra",
#   controls = c("police", "income", "poverty", "unemployrt"),
#   aggregate = "cohort"
# )

## ----vce----------------------------------------------------------------------
# HC1 is recommended for the smoking data (single treated unit makes HC3
# numerically unstable due to hat value = 1)
result_hc1 <- suppressWarnings(lwdid(
  data = smoking,
  y = "lcigsale",
  ivar = "state",
  tvar = "year",
  d = "d",
  post = "post",
  rolling = "demean",
  estimator = "ra",
  vce = "hc1"
))

summary(result_hc1)

## ----extract------------------------------------------------------------------
coef(result)
confint(result)
vcov(result)
nobs(result)

## ----parallel, eval=FALSE-----------------------------------------------------
# result_par <- lwdid(
#   data = castle,
#   y = "lhomicide",
#   ivar = "sid",
#   tvar = "year",
#   gvar = "gvar",
#   rolling = "demean",
#   estimator = "ra",
#   aggregate = "cohort",
#   parallel = TRUE,
#   n_cores = 4
# )

