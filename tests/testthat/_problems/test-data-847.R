# Extracted from test-data.R:847

# prequel ----------------------------------------------------------------------
load("../../data/smoking.rda")
load("../../data/castle.rda")
source("../../R/simulate.R")

# test -------------------------------------------------------------------------
skip_if_not(exists("lwdid", mode = "function"),
              message = "lwdid() not yet available (Epic 2+)")
res <- lwdid(data = castle, yname = "lhomicide", idname = "sid",
               tname = "year", gname = "gvar",
               xformla = ~ cdl, est_method = "ra", vce = "hc1",
               transform = "detrend")
