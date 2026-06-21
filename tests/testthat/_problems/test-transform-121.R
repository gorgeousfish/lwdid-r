# Extracted from test-transform.R:121

# test -------------------------------------------------------------------------
dt <- data.table::data.table(
    id = c(rep(1, 6), rep(2, 6), 3, 3, 3, 3),
    t  = c(1:6, 1:6, 3, 4, 5, 6),
    y  = c(1, 2, 3, 4, 5, 6,
           2, 3, 4, 5, 6, 7,
           10, 11, 12, 13)
  )
result <- expect_warning(
    suppressMessages(
      transform_detrend(dt, "y", "id", "t", g = 4L)
    ),
    class = "lwdid_data"
  )
u3 <- result[result[["id"]] == 3]
expect_equal(u3[["slope"]][1], 0)
u3_t4 <- result[result[["id"]] == 3 & result[["t"]] == 4]
u3_t5 <- result[result[["id"]] == 3 & result[["t"]] == 5]
u3_t6 <- result[result[["id"]] == 3 & result[["t"]] == 6]
expect_equal(u3_t4[["y_trans"]], 11 - 10, tolerance = 1e-15)
expect_equal(u3_t5[["y_trans"]], 12 - 10, tolerance = 1e-15)
expect_equal(u3_t6[["y_trans"]], 13 - 10, tolerance = 1e-15)
u1 <- result[result[["id"]] == 1]
u2 <- result[result[["id"]] == 2]
expect_true(u1[["slope"]][1] != 0)
