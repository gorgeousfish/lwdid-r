# Extracted from test-transform.R:180

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
