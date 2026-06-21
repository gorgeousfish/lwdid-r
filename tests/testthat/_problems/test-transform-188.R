# Extracted from test-transform.R:188

# test -------------------------------------------------------------------------
dt <- data.table::data.table(
    id = c(1, 1, 1, 1, 1, 2, 2, 2, 2),
    t  = c(1, 2, 3, 4, 5, 3, 3, 4, 5),
    y  = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  )
result <- expect_warning(
    suppressMessages(
      transform_detrend(dt, "y", "id", "t", g = 4L)
    ),
    class = "lwdid_data"
  )
u2 <- result[result[["id"]] == 2]
expect_equal(u2[["slope"]][1], 0)
pre_mean_u2 <- mean(c(6, 7))
u2_t4 <- result[result[["id"]] == 2 & result[["t"]] == 4]
expect_equal(
    u2_t4[["y_trans"]][1], 8 - pre_mean_u2, tolerance = 1e-15
  )
