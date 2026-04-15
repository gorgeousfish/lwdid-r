# ============================================================================
# test-data.R — Tests for built-in datasets and simulate_panel_data()
#
# Test Groups:
#   A: smoking basic validation
#   B: smoking panel structure
#   C: smoking numerical range
#   D: castle basic validation
#   E: castle panel structure
#   F: castle cohort structure
#   G: castle numerical range
#   G+: castle cohort state name verification
#   G-ATT: castle ATT benchmarks (conditional skip)
#   H: Python/Stata consistency (conditional skip)
#   I: simulate_panel_data() tests
# ============================================================================

# ==========================================================================
# Group A — smoking basic validation
# ==========================================================================

test_that("A1: smoking loads as data.frame", {
  data(smoking, package = "lwdid")
  expect_true(is.data.frame(smoking))
})

test_that("A2: smoking has correct dimensions", {
  data(smoking, package = "lwdid")
  expect_equal(nrow(smoking), 1209L)
  expect_equal(ncol(smoking), 12L)
})

test_that("A3: smoking has exact column names", {
  data(smoking, package = "lwdid")
  expected_cols <- c(
    "state", "year", "cigsale", "lcigsale",
    "retprice", "lretprice", "lnincome", "beer",
    "age15to24", "d", "post", "treat"
  )
  expect_equal(names(smoking), expected_cols)
})

test_that("A4: smoking column types are correct", {
  data(smoking, package = "lwdid")
  expect_true(is.integer(smoking$state))
  expect_true(is.integer(smoking$year))
  expect_true(is.numeric(smoking$cigsale))
  expect_true(is.numeric(smoking$lcigsale))
  expect_true(is.numeric(smoking$retprice))
  expect_true(is.numeric(smoking$lretprice))
  expect_true(is.numeric(smoking$lnincome))
  expect_true(is.numeric(smoking$beer))
  expect_true(is.numeric(smoking$age15to24))
  expect_true(is.integer(smoking$d))
  expect_true(is.integer(smoking$post))
  expect_true(is.integer(smoking$treat))
})

test_that("A5: smoking key columns have no NA", {
  data(smoking, package = "lwdid")
  expect_false(anyNA(smoking$d))
  expect_false(anyNA(smoking$post))
  expect_false(anyNA(smoking$lcigsale))
  expect_false(anyNA(smoking$state))
  expect_false(anyNA(smoking$year))
})

test_that("A6: smoking allows NA in optional columns", {
  data(smoking, package = "lwdid")
  expect_true("lnincome" %in% names(smoking))
  expect_true("beer" %in% names(smoking))
  expect_true("age15to24" %in% names(smoking))
})

# ==========================================================================
# Group B — smoking panel structure
# ==========================================================================

test_that("B1: smoking has 39 unique states", {
  data(smoking, package = "lwdid")
  expect_equal(length(unique(smoking$state)), 39L)
})

test_that("B2: smoking has 31 unique years (1970-2000)", {
  data(smoking, package = "lwdid")
  expect_equal(length(unique(smoking$year)), 31L)
  expect_equal(range(smoking$year), c(1970L, 2000L))
})

test_that("B3: smoking is a balanced panel", {
  data(smoking, package = "lwdid")
  obs_per_state <- table(smoking$state)
  expect_true(all(obs_per_state == 31L))
})

test_that("B4: no duplicate (state, year) pairs", {
  data(smoking, package = "lwdid")
  pairs <- paste(smoking$state, smoking$year, sep = "-")
  expect_equal(length(pairs), length(unique(pairs)))
})

test_that("B5: only California (state=3) is treated", {
  data(smoking, package = "lwdid")
  treated_states <- unique(smoking$state[smoking$d == 1L])
  expect_equal(treated_states, 3L)
  expect_equal(sum(smoking$d == 1L), 31L)
})

test_that("B6: 38 control states have d=0", {
  data(smoking, package = "lwdid")
  control_states <- unique(smoking$state[smoking$d == 0L])
  expect_equal(length(control_states), 38L)
})

test_that("B7: post=1 for years 1989-2000", {
  data(smoking, package = "lwdid")
  post_years <- sort(unique(smoking$year[smoking$post == 1L]))
  expect_equal(post_years, 1989L:2000L)
})

test_that("B8: treat = d * post (element-wise)", {
  data(smoking, package = "lwdid")
  expect_identical(smoking$treat, smoking$d * smoking$post)
})

# ==========================================================================
# Group C — smoking numerical range
# ==========================================================================

test_that("C1: smoking state and year ranges", {
  data(smoking, package = "lwdid")
  expect_true(all(smoking$state >= 1L & smoking$state <= 39L))
  expect_true(all(smoking$year >= 1970L & smoking$year <= 2000L))
})

test_that("C2: smoking lcigsale in reasonable range", {
  data(smoking, package = "lwdid")
  expect_true(all(smoking$lcigsale >= 3.5))
  expect_true(all(smoking$lcigsale <= 5.8))
})

test_that("C3: smoking d and post are binary", {
  data(smoking, package = "lwdid")
  expect_true(all(smoking$d %in% c(0L, 1L)))
  expect_true(all(smoking$post %in% c(0L, 1L)))
})

test_that("C4: lcigsale approx log(cigsale)", {
  data(smoking, package = "lwdid")
  expect_equal(
    smoking$lcigsale,
    log(smoking$cigsale),
    tolerance = 1e-4
  )
})

test_that("C5: PRD 14 column is d not treated", {
  data(smoking, package = "lwdid")
  expect_true("d" %in% names(smoking))
  expect_false("treated" %in% names(smoking))
})

test_that("C6: treat column exists and equals d * post", {
  data(smoking, package = "lwdid")
  expect_true("treat" %in% names(smoking))
  expect_identical(smoking$treat, smoking$d * smoking$post)
})
# Group D — castle basic
test_that("D1: castle is data.frame", {
  expect_true(is.data.frame(castle))
})
test_that("D2: castle dims 550x18", {
  expect_equal(nrow(castle), 550L)
  expect_equal(ncol(castle), 18L)
})
test_that("D3: castle colnames", {
  ec <- c("sid","state_name","year","lhomicide","homicide",
    "population","police","unemployrt","income",
    "poverty","blackm_15_24","whitem_15_24","prisoner",
    "cdl","post","effyear","gvar","dinf")
  expect_identical(names(castle), ec)
})
test_that("D4: castle coltypes", {
  expect_type(castle$sid, "integer")
  expect_type(castle$year, "integer")
  expect_type(castle$post, "integer")
  expect_type(castle$effyear, "integer")
  expect_type(castle$gvar, "integer")
  expect_type(castle$dinf, "integer")
  expect_type(castle$state_name, "character")
  expect_type(castle$lhomicide, "double")
  expect_type(castle$homicide, "double")
  expect_type(castle$population, "double")
  expect_type(castle$police, "double")
  expect_type(castle$unemployrt, "double")
  expect_type(castle$income, "double")
  expect_type(castle$poverty, "double")
  expect_type(castle$blackm_15_24, "double")
  expect_type(castle$whitem_15_24, "double")
  expect_type(castle$prisoner, "double")
  expect_type(castle$cdl, "double")
})
test_that("D5: lhomicide no NA", {
  expect_false(anyNA(castle$lhomicide))
})
# Group E — castle panel structure
test_that("E1: 50 unique states", {
  expect_equal(length(unique(castle$sid)), 50L)
})
test_that("E2: 11 years 2000-2010", {
  uy <- sort(unique(castle$year))
  expect_equal(length(uy), 11L)
  expect_identical(uy, 2000L:2010L)
})
test_that("E3: balanced 11 obs/state", {
  expect_true(all(table(castle$sid) == 11L))
})
test_that("E4: no dup (sid,year)", {
  k <- paste(castle$sid, castle$year, sep="_")
  expect_equal(length(k), length(unique(k)))
})
test_that("E5: sid skips 9", {
  expect_false(9L %in% castle$sid)
})
# Group F — castle cohort structure
test_that("F1: 29 NT states (gvar=NA)", {
  gpu <- tapply(castle$gvar, castle$sid, function(x) x[1])
  expect_equal(sum(is.na(gpu)), 29L)
})
test_that("F2: 21 treated states", {
  gpu <- tapply(castle$gvar, castle$sid, function(x) x[1])
  expect_equal(sum(!is.na(gpu)), 21L)
})
test_that("F3: cohort sizes exact", {
  gpu <- tapply(castle$gvar, castle$sid, function(x) x[1])
  tg <- gpu[!is.na(gpu)]
  cc <- table(tg)
  expect_equal(as.integer(cc["2005"]), 1L)
  expect_equal(as.integer(cc["2006"]), 13L)
  expect_equal(as.integer(cc["2007"]), 4L)
  expect_equal(as.integer(cc["2008"]), 2L)
  expect_equal(as.integer(cc["2009"]), 1L)
  expect_equal(length(cc), 5L)
})
test_that("F4: gvar==effyear consistent", {
  nn <- !is.na(castle$gvar)
  expect_true(all(castle$gvar[nn] == castle$effyear[nn]))
  expect_identical(is.na(castle$gvar), is.na(castle$effyear))
})
test_that("F5: dinf==1 iff gvar NA", {
  expect_true(all((castle$dinf==1L) == is.na(castle$gvar)))
})
test_that("F6: gvar non-NA in {2005..2009}", {
  nng <- castle$gvar[!is.na(castle$gvar)]
  expect_true(all(nng %in% 2005L:2009L))
})
# Group G — castle numerical range
test_that("G1: year [2000,2010]", {
  expect_equal(min(castle$year), 2000L)
  expect_equal(max(castle$year), 2010L)
})
test_that("G2: lhomicide in [-1,3]", {
  expect_true(min(castle$lhomicide) >= -1.0)
  expect_true(max(castle$lhomicide) <= 3.0)
})
test_that("G3: cdl in [0,1]", {
  expect_true(min(castle$cdl) >= 0.0)
  expect_true(max(castle$cdl) <= 1.0)
})
test_that("G4: cdl continuous", {
  expect_true(any(castle$cdl > 0 & castle$cdl < 1))
})
test_that("G5: lhomicide approx log(homicide)", {
  expect_equal(as.numeric(castle$lhomicide),
    log(as.numeric(castle$homicide)), tolerance=1e-4)
})
# Group G+ — cohort state names
test_that("G+1: 2005=Florida", {
  s <- unique(castle$state_name[castle$gvar==2005L & !is.na(castle$gvar)])
  expect_true("Florida" %in% s)
})
test_that("G+2: 2007=MO,ND,TN,TX", {
  s <- unique(castle$state_name[castle$gvar==2007L & !is.na(castle$gvar)])
  expect_true("Missouri" %in% s)
  expect_true("North Dakota" %in% s)
  expect_true("Tennessee" %in% s)
  expect_true("Texas" %in% s)
})
test_that("G+3: 2008=OH,WV", {
  s <- unique(castle$state_name[castle$gvar==2008L & !is.na(castle$gvar)])
  expect_true("Ohio" %in% s)
  expect_true("West Virginia" %in% s)
})
test_that("G+4: 2009=Montana", {
  s <- unique(castle$state_name[castle$gvar==2009L & !is.na(castle$gvar)])
  expect_true("Montana" %in% s)
})
# Group G-ATT — castle ATT benchmarks (conditional)
test_that("G-ATT1: demean ATT~0.092", {
  skip("ATT benchmarks require full estimation engine (Epic 2+)")
  data(castle, package = "lwdid")
  res <- lwdid(data=castle, yname="lhomicide", idname="sid",
    tname="year", gname="gvar",
    xformla=~cdl, est_method="ra", vce="hc1")
  expect_equal(res$att, 0.092, tolerance=0.03)
})
test_that("G-ATT2: detrend ATT~0.067", {
  skip("ATT benchmarks require full estimation engine (Epic 2+)")
  data(castle, package = "lwdid")
  res <- lwdid(data=castle, yname="lhomicide", idname="sid",
    tname="year", gname="gvar",
    xformla=~cdl, est_method="ra", vce="hc1",
    transform="detrend")
  expect_equal(res$att, 0.067, tolerance=0.03)
})
# Group H — Python/Stata consistency (optional)
test_that("H1: smoking vs Python CSV", {
  skip("External reference data not available")
  py <- read.csv("../../../lwdid-py_v0.2.3/data/smoking.csv")
  expect_equal(nrow(py), 1209L)
})
test_that("H2: castle vs Python CSV", {
  skip("External reference data not available")
  py <- read.csv("../../../lwdid-py_v0.2.3/data/castle.csv")
  expect_equal(nrow(py), 550L)
})
# Group I — simulate_panel_data() tests
# I.1 Common Timing basic
test_that("I1a: CT columns", {
  ct <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=42)
  expect_true(is.data.frame(ct))
  expect_identical(names(ct), c("id","year","y","d","post"))
})
test_that("I1b: CT dimensions", {
  ct <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=42)
  expect_equal(nrow(ct), 100L)
  expect_equal(ncol(ct), 5L)
})
test_that("I1c: CT d split correct", {
  ct <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=42)
  dpu <- tapply(ct$d, ct$id, function(x) x[1])
  expect_equal(sum(dpu==1L), 8L)
  expect_equal(sum(dpu==0L), 12L)
})
test_that("I1d: CT post switches at correct year", {
  ct <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=42)
  pre <- sort(unique(ct$year[ct$post==0L]))
  post <- sort(unique(ct$year[ct$post==1L]))
  expect_identical(pre, c(2001L, 2002L))
  expect_identical(post, c(2003L, 2004L, 2005L))
})
test_that("I1e: CT years 2001..2000+n_periods", {
  ct <- simulate_panel_data(n_units=10, n_periods=8,
    n_treated=3, treatment_period=4, seed=1)
  expect_identical(sort(unique(ct$year)), 2001L:2008L)
})
# I.2 Staggered mode basic
test_that("I2a: Stag columns", {
  sg <- simulate_panel_data(n_units=60, n_periods=10,
    cohorts=list("2005"=20,"2007"=20),
    n_never_treated=20, seed=42)
  expect_true(is.data.frame(sg))
  expect_identical(names(sg), c("id","year","y","gvar"))
})
test_that("I2b: Stag dimensions", {
  sg <- simulate_panel_data(n_units=60, n_periods=10,
    cohorts=list("2005"=20,"2007"=20),
    n_never_treated=20, seed=42)
  expect_equal(nrow(sg), 600L)
  expect_equal(ncol(sg), 4L)
})
test_that("I2c: Stag cohort structure correct", {
  sg <- simulate_panel_data(n_units=60, n_periods=10,
    cohorts=list("2005"=20,"2007"=20),
    n_never_treated=20, seed=42)
  gpu <- tapply(sg$gvar, sg$id, function(x) x[1])
  expect_equal(sum(gpu==2005L, na.rm=TRUE), 20L)
  expect_equal(sum(gpu==2007L, na.rm=TRUE), 20L)
  expect_equal(sum(is.na(gpu)), 20L)
})
test_that("I2d: Stag NT units have gvar=NA (not 0)", {
  sg <- simulate_panel_data(n_units=60, n_periods=10,
    cohorts=list("2005"=20,"2007"=20),
    n_never_treated=20, seed=42)
  nt_rows <- sg[sg$id > 40, ]
  expect_true(all(is.na(nt_rows$gvar)))
  expect_false(any(nt_rows$gvar == 0L, na.rm=TRUE))
})
# I.3 Seed reproducibility
test_that("I3a: same seed same data", {
  d1 <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=123)
  d2 <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=123)
  expect_identical(d1, d2)
})
test_that("I3b: different seed different data", {
  d1 <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=123)
  d2 <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=456)
  expect_false(identical(d1$y, d2$y))
})
# I.4 Function-valued ATT
test_that("I4a: function ATT produces varying effects", {
  fn <- function(i, t, g) 0.5 * (t - g + 1)
  ct <- simulate_panel_data(n_units=30, n_periods=8,
    n_treated=10, treatment_period=4, true_att=fn, seed=42)
  tr_post <- ct[ct$d==1L & ct$post==1L, ]
  y_by_year <- tapply(tr_post$y, tr_post$year, mean)
  diffs <- diff(as.numeric(y_by_year))
  expect_true(all(diffs > 0))
})
test_that("I4b: function ATT in staggered mode", {
  fn <- function(i, t, g) 0.5 * (t - g + 1)
  sg <- simulate_panel_data(n_units=40, n_periods=10,
    cohorts=list("2004"=15,"2007"=15),
    n_never_treated=10, true_att=fn, seed=42)
  expect_equal(nrow(sg), 400L)
  expect_true(is.data.frame(sg))
})
# I.5 Unbalanced panel
test_that("I5a: unbalanced has fewer rows", {
  bal <- simulate_panel_data(n_units=30, n_periods=8,
    n_treated=10, treatment_period=4, seed=42)
  unb <- simulate_panel_data(n_units=30, n_periods=8,
    n_treated=10, treatment_period=4,
    balanced=FALSE, seed=42)
  expect_true(nrow(unb) < nrow(bal))
})
test_that("I5b: unbalanced min 2 obs per unit", {
  unb <- simulate_panel_data(n_units=30, n_periods=8,
    n_treated=10, treatment_period=4,
    balanced=FALSE, seed=42)
  obs <- table(unb$id)
  expect_true(all(obs >= 2L))
})
test_that("I5c: unbalanced varying obs counts", {
  unb <- simulate_panel_data(n_units=30, n_periods=8,
    n_treated=10, treatment_period=4,
    balanced=FALSE, seed=42)
  obs <- table(unb$id)
  expect_true(length(unique(as.integer(obs))) > 1L)
})
# I.6 Unit-specific trend
test_that("I6: trend parameter adds variation", {
  set.seed(99)
  d_no <- simulate_panel_data(n_units=50, n_periods=10,
    n_treated=20, treatment_period=5,
    trend=NULL, seed=99)
  d_tr <- simulate_panel_data(n_units=50, n_periods=10,
    n_treated=20, treatment_period=5,
    trend=0.3, seed=99)
  var_no <- var(d_no$y)
  var_tr <- var(d_tr$y)
  expect_true(var_tr > var_no)
})
# I.7 Edge cases
test_that("I7a: n_units=3 works", {
  d <- simulate_panel_data(n_units=3, n_periods=4,
    n_treated=1, treatment_period=2, seed=1)
  expect_equal(nrow(d), 12L)
  expect_equal(ncol(d), 5L)
})
test_that("I7b: n_controls=0 no x columns", {
  d <- simulate_panel_data(n_units=10, n_periods=5,
    n_treated=3, treatment_period=3,
    n_controls=0, seed=1)
  expect_identical(names(d), c("id","year","y","d","post"))
})
test_that("I7c: n_controls=2 adds x1,x2", {
  d <- simulate_panel_data(n_units=10, n_periods=5,
    n_treated=3, treatment_period=3,
    n_controls=2, seed=1)
  expect_true("x1" %in% names(d))
  expect_true("x2" %in% names(d))
  expect_equal(ncol(d), 7L)
})
# I.8 Mode dispatch errors
test_that("I8a: no args errors", {
  expect_error(simulate_panel_data(),
    "Must specify either")
})
test_that("I8b: both modes errors", {
  expect_error(
    simulate_panel_data(n_treated=10,
      treatment_period=3,
      cohorts=list("2005"=10)),
    "Cannot specify both")
})
# I.9 Cohorts year-key convention
test_that("I9a: year keys produce correct gvar", {
  sg <- simulate_panel_data(n_units=50, n_periods=10,
    cohorts=list("2006"=20,"2008"=20),
    n_never_treated=10, seed=42)
  gpu <- tapply(sg$gvar, sg$id, function(x) x[1])
  expect_equal(sum(gpu==2006L, na.rm=TRUE), 20L)
  expect_equal(sum(gpu==2008L, na.rm=TRUE), 20L)
})
test_that("I9b: D_it switches at correct year", {
  sg <- simulate_panel_data(n_units=30, n_periods=10,
    cohorts=list("2006"=10,"2008"=10),
    n_never_treated=10, true_att=5.0,
    unit_fe_sd=0, time_trend_coef=0,
    error_sd=0, seed=42)
  u1 <- sg[sg$id==1, ]
  expect_true(u1$gvar[1] == 2006L)
  pre <- u1$y[u1$year < 2006L]
  post <- u1$y[u1$year >= 2006L]
  expect_true(all(abs(pre - 0) < 1e-10))
  expect_true(all(abs(post - 5.0) < 1e-10))
})
# I.10 Large sample ATT recovery
test_that("I10: large sample recovers ATT", {
  d <- simulate_panel_data(n_units=1000, n_periods=10,
    n_treated=400, treatment_period=6,
    true_att=2.0, unit_fe_sd=1.0,
    error_sd=1.0, seed=42)
  tr_post <- d[d$d==1L & d$post==1L, ]
  tr_pre <- d[d$d==1L & d$post==0L, ]
  co_post <- d[d$d==0L & d$post==1L, ]
  co_pre <- d[d$d==0L & d$post==0L, ]
  did_est <- (mean(tr_post$y) - mean(tr_pre$y)) -
    (mean(co_post$y) - mean(co_pre$y))
  expect_equal(did_est, 2.0, tolerance=0.3)
})
# I.11 lwdid() compatibility
test_that("I11: CT output compatible with lwdid", {
  ct <- simulate_panel_data(n_units=20, n_periods=5,
    n_treated=8, treatment_period=3, seed=42)
  expect_true(is.integer(ct$id))
  expect_true(is.integer(ct$year))
  expect_true(is.numeric(ct$y))
  expect_true(is.integer(ct$d))
  expect_true(is.integer(ct$post))
})
