# Extracted from test-transform.R:157

# test -------------------------------------------------------------------------
dt <- data.table::data.table(
    id = rep(1, 4),
    t  = 1:4,
    y  = c(2, 4, 7, 9)
  )
result <- expect_warning(
    suppressMessages(
      transform_detrend(dt, "y", "id", "t", g = 3L)
    ),
    class = "lwdid_small_sample"
  )
u1 <- result[result[["id"]] == 1]
expect_equal(u1[["intercept_c"]][1], 3, tolerance = 1e-10)
expect_equal(u1[["slope"]][1], 2, tolerance = 1e-10)
expect_equal(u1[["t_bar_pre"]][1], 1.5, tolerance = 1e-10)
u1_t3 <- result[result[["id"]] == 1 & result[["t"]] == 3]
u1_t4 <- result[result[["id"]] == 1 & result[["t"]] == 4]
expect_equal(u1_t3[["y_trans"]], 1, tolerance = 1e-10)
