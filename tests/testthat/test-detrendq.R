# test-detrendq.R
# Story E9-03 dedicated detrendq regressions

make_e903_detrendq_panel <- function() {
  seasonal_effects <- c(`1` = 0, `2` = 1, `3` = 0.5, `4` = 1.5)
  quarter_pattern <- c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L)
  post_pattern <- c(rep.int(0L, 6L), rep.int(1L, 4L))
  intercepts <- c(5, 15, 25, 35)
  slopes <- c(0.5, 0.7, 0.3, 0.4)

  do.call(
    rbind,
    lapply(seq_along(intercepts), function(idx) {
      effect <- if (idx <= 2L) 4 else 0
      data.frame(
        id = idx,
        time = 1:10,
        quarter = quarter_pattern,
        d = as.integer(idx <= 2L),
        post = post_pattern,
        y = intercepts[[idx]] +
          slopes[[idx]] * (1:10) +
          unname(seasonal_effects[as.character(quarter_pattern)]) +
          effect * post_pattern
      )
    })
  )
}

make_e903_detrendq_exclude_surface_panel <- function() {
  seasonal_effects <- c(`1` = 0, `2` = 1, `3` = 0.5, `4` = 1.5)
  quarter_pattern <- c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L)
  post_pattern <- c(rep.int(0L, 7L), rep.int(1L, 3L))
  intercepts <- c(5, 15, 25, 35)
  slopes <- c(0.5, 0.7, 0.3, 0.4)
  anticipation_bumps <- list(
    c(0, 0, 0, 0, 0, 0, 2, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, -1.5, 0, 0, 0),
    rep.int(0, 10),
    rep.int(0, 10)
  )

  do.call(
    rbind,
    lapply(seq_along(intercepts), function(idx) {
      effect <- if (idx <= 2L) 4 else 0
      data.frame(
        id = idx,
        time = 1:10,
        quarter = quarter_pattern,
        d = as.integer(idx <= 2L),
        post = post_pattern,
        y = intercepts[[idx]] +
          slopes[[idx]] * (1:10) +
          unname(seasonal_effects[as.character(quarter_pattern)]) +
          effect * post_pattern +
          anticipation_bumps[[idx]]
      )
    })
  )
}

test_that("TC-9.3.1: transform_detrendq() reproduces the centered trend-plus-season fit", {
  data <- make_e903_detrendq_panel()

  result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )

  expect_equal(
    result[result$id == 1L, "seasonal_fit"][[1]],
    c(5.5, 7, 7, 8.5, 7.5, 9, 9, 10.5, 9.5, 11),
    tolerance = 1e-10
  )
  expect_equal(
    result[result$id == 1L, "y_trans"][[1]],
    c(rep.int(0, 6L), rep.int(4, 4L)),
    tolerance = 1e-10
  )
  expect_identical(unique(result[result$id == 1L, "n_pre"][[1]]), 6L)
  expect_equal(unique(result[result$id == 1L, "intercept_c"][[1]]), 6.75, tolerance = 1e-10)
  expect_equal(unique(result[result$id == 1L, "slope"][[1]]), 0.5, tolerance = 1e-10)
  expect_equal(unique(result[result$id == 1L, "t_bar_pre"][[1]]), 3.5, tolerance = 1e-10)
})

test_that("TC-9.3.2: transform_common() dispatches to detrendq and preserves common-timing residualization", {
  data <- make_e903_detrendq_panel()

  result <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )

  expect_equal(
    result[result$id == 1L, "y_trans"][[1]],
    c(rep.int(0, 6L), rep.int(4, 4L)),
    tolerance = 1e-10
  )
  expect_equal(
    result[result$id == 3L, "y_trans"][[1]],
    rep.int(0, 10L),
    tolerance = 1e-10
  )
})

test_that("TC-9.3.3: top-level lwdid detrendq returns the expected common-timing public surface", {
  data <- make_e903_detrendq_panel()

  result <- suppressWarnings(
    suppressMessages(
      lwdid(
        data = data,
        y = "y",
        ivar = "id",
        tvar = "time",
        d = "d",
        post = "post",
        rolling = "detrendq",
        season_var = "quarter",
        Q = 4L
      )
    )
  )

  expect_s3_class(result, "lwdid_result")
  expect_identical(result$rolling, "detrendq")
  expect_equal(result$att, 4, tolerance = 1e-10)
  expect_identical(result$K, 6L)
  expect_identical(result$tpost1, 7L)
})

test_that("TC-9.3.4: detrendq quarter alias matches season_var in transform_common()", {
  data <- make_e903_detrendq_panel()

  baseline <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )
  aliased <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    rolling = "detrendq",
    post = "post",
    quarter = "quarter",
    Q = 4L
  )

  expect_equal(aliased$y_trans, baseline$y_trans, tolerance = 1e-10)
  expect_equal(aliased$seasonal_fit, baseline$seasonal_fit, tolerance = 1e-10)
  expect_equal(aliased$n_pre, baseline$n_pre, tolerance = 0)
})

test_that("TC-9.3.5: detrendq exclude_pre_periods changes the dedicated seasonal surface and top-level setting", {
  data <- make_e903_detrendq_exclude_surface_panel()

  baseline <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 8L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L,
    exclude_pre_periods = 0L
  )
  excluded <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 8L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L,
    exclude_pre_periods = 1L
  )

  expect_identical(unique(baseline[baseline$id == 1L, "n_pre"][[1]]), 7L)
  expect_identical(unique(excluded[excluded$id == 1L, "n_pre"][[1]]), 6L)
  expect_false(isTRUE(all.equal(
    baseline[baseline$id == 1L, "seasonal_fit"][[1]],
    excluded[excluded$id == 1L, "seasonal_fit"][[1]],
    tolerance = 1e-10
  )))
  expect_false(isTRUE(all.equal(
    baseline[baseline$id == 1L, "y_trans"][[1]],
    excluded[excluded$id == 1L, "y_trans"][[1]],
    tolerance = 1e-10
  )))

  baseline_result <- suppressWarnings(
    suppressMessages(
      lwdid(
        data = data,
        y = "y",
        ivar = "id",
        tvar = "time",
        d = "d",
        post = "post",
        rolling = "detrendq",
        season_var = "quarter",
        Q = 4L,
        exclude_pre_periods = 0L
      )
    )
  )
  excluded_result <- suppressWarnings(
    suppressMessages(
      lwdid(
        data = data,
        y = "y",
        ivar = "id",
        tvar = "time",
        d = "d",
        post = "post",
        rolling = "detrendq",
        season_var = "quarter",
        Q = 4L,
        exclude_pre_periods = 1L
      )
    )
  )

  expect_identical(baseline_result$exclude_pre_periods, 0L)
  expect_identical(excluded_result$exclude_pre_periods, 1L)
  expect_true(is.finite(baseline_result$att))
  expect_true(is.finite(excluded_result$att))
  expect_false(isTRUE(all.equal(
    baseline_result$att,
    excluded_result$att,
    tolerance = 1e-10
  )))
})
