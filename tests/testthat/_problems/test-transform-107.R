# Extracted from test-transform.R:107

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
