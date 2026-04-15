# test-lwdid-estimators.R — E2E integration tests for lwdid() estimators
# Task E6-05.5
if (!exists(".lwdid_env") || is.null(.lwdid_env$warning_registry)) {
  .lwdid_env$warning_registry <- new_warning_registry()
}
# CT data: 8 periods, d/post columns for common_timing path, true ATT ~ 2.0
make_ct_data <- function(n = 300, seed = 123) {
  set.seed(seed)
  n_treat <- floor(n / 2)
  ids <- seq_len(n)
  years <- 2001:2008
  df <- expand.grid(id = ids, year = years)
  df <- df[order(df$id, df$year), ]
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  df$x1 <- x1[df$id]
  df$x2 <- x2[df$id]
  df$d <- as.integer(df$id <= n_treat)
  df$post <- as.integer(df$year > 2004)
  df$y <- 1 + 0.5 * df$x1 - 0.3 * df$x2 + 0.3 * (df$year - 2001) +
    2.0 * (df$d == 1 & df$post == 1) + rnorm(nrow(df), sd = 0.5)
  df
}
# Staggered data: 10 periods, gvar column, true ATT ~ 2.0
make_stag_data <- function(n = 100, seed = 456) {
  set.seed(seed)
  n_nt <- floor(n / 3)
  n_c1 <- floor(n / 3)
  n_c2 <- n - n_nt - n_c1
  ids <- seq_len(n)
  years <- 2001:2010
  df <- expand.grid(id = ids, year = years)
  df <- df[order(df$id, df$year), ]
  gvar_map <- c(rep(0L, n_nt), rep(2005L, n_c1), rep(2007L, n_c2))
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  df$gvar <- gvar_map[df$id]
  df$x1 <- x1[df$id]
  df$x2 <- x2[df$id]
  treated <- as.integer(df$gvar > 0 & df$year >= df$gvar)
  df$y <- 1 + 0.5 * df$x1 - 0.3 * df$x2 + 0.3 * (df$year - 2001) +
    2.0 * treated + rnorm(nrow(df), sd = 0.5)
  df
}
check_ct <- function(r, lab = "") {
  expect_true(is.numeric(r$att) && is.finite(r$att), info = paste(lab, "att"))
  expect_true(r$se_att > 0, info = paste(lab, "se"))
  expect_true(is.finite(r$ci_lower), info = paste(lab, "ci_lo"))
  expect_true(is.finite(r$ci_upper), info = paste(lab, "ci_hi"))
  expect_true(r$ci_lower < r$ci_upper, info = paste(lab, "ci_ord"))
}
