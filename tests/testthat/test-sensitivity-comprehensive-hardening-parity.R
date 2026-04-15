library(testthat)

resolve_hardening_audit_path <- function(filename) {
  resolve_parity_fixture_path(filename)
}

resolve_smoking_raw_csv_for_hardening <- function() {
  resolve_parity_fixture_path("smoking.csv")
}

load_smoking_frozen_controls_for_hardening <- function() {
  raw_csv_path <- resolve_smoking_raw_csv_for_hardening()
  fixture_path <- resolve_hardening_audit_path(
    "e8_03_smoking_frozen_controls_fixture.csv"
  )

  smoking_raw <- utils::read.csv(raw_csv_path, stringsAsFactors = FALSE)
  fixture <- utils::read.csv(fixture_path, stringsAsFactors = FALSE)
  controls <- c("lnincome", "beer", "age15to24", "lretprice")

  smoking_frozen <- merge(
    smoking_raw[, setdiff(names(smoking_raw), controls), drop = FALSE],
    fixture,
    by = "state",
    all.x = TRUE,
    sort = FALSE
  )

  smoking_frozen <- smoking_frozen[
    order(match(smoking_frozen$state, unique(smoking_raw$state)), smoking_frozen$year),
    ,
    drop = FALSE
  ]

  list(
    data = smoking_frozen,
    controls = controls
  )
}

collect_numeric_scalars <- function(x, prefix = NULL) {
  out <- list()

  if (is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      nms <- rep("", length(x))
    }

    for (i in seq_along(x)) {
      child_name <- nms[[i]]
      if (!nzchar(child_name)) {
        child_name <- sprintf("[[%d]]", i)
      }
      child_prefix <- if (is.null(prefix)) {
        child_name
      } else {
        paste0(prefix, "$", child_name)
      }
      out <- c(out, collect_numeric_scalars(x[[i]], child_prefix))
    }

    return(out)
  }

  if (is.numeric(x) && length(x) == 1L) {
    out[[prefix]] <- unname(as.numeric(x))
  }

  out
}

is_att_like_path <- function(path) {
  grepl("\\$att$", path) ||
    grepl(
      "^transformation_comparison\\$(demean_att|detrend_att)$",
      path
    ) ||
    grepl(
      "^estimator_comparison\\$(ra|ipw|ipwra|baseline_att)$",
      path
    )
}

compute_numeric_audit <- function(result, y_sd, att_bound_multiplier = 10) {
  numeric_map <- collect_numeric_scalars(result)
  numeric_values <- unlist(numeric_map, use.names = TRUE)

  nonfinite_paths <- names(numeric_values)[!is.finite(numeric_values)]
  att_paths <- names(numeric_values)[vapply(names(numeric_values), is_att_like_path, logical(1))]
  att_values <- numeric_values[att_paths]
  att_ratio_by_path <- abs(att_values) / y_sd
  out_of_range_paths <- names(att_ratio_by_path)[att_ratio_by_path >= att_bound_multiplier]

  list(
    numeric_field_count = length(numeric_values),
    nonfinite_paths = sort(unname(nonfinite_paths)),
    att_ratio_by_path = att_ratio_by_path,
    max_att_ratio = if (length(att_ratio_by_path) == 0L) {
      0
    } else {
      max(att_ratio_by_path)
    },
    out_of_range_paths = sort(unname(out_of_range_paths))
  )
}

normalize_oracle_character <- function(x) {
  if (length(x) == 0L) {
    return(character(0))
  }

  sort(unlist(x, use.names = FALSE))
}

normalize_oracle_att_ratios <- function(x) {
  if (length(x) == 0L) {
    return(setNames(numeric(0), character(0)))
  }

  unlist(x, use.names = TRUE)
}

run_hardening_case <- function(case_id) {
  warnings_seen <- character(0)

  case <- switch(
    case_id,
    "smoking-no-controls" = {
      data("smoking", package = "lwdid", envir = environment())
      list(
        data = smoking,
        args = list(
          y = "lcigsale",
          ivar = "state",
          tvar = "year",
          d = "d",
          post = "post",
          controls = NULL
        )
      )
    },
    "castle-no-controls" = {
      data("castle", package = "lwdid", envir = environment())
      list(
        data = castle,
        args = list(
          y = "lhomicide",
          ivar = "sid",
          tvar = "year",
          gvar = "gvar",
          controls = NULL
        )
      )
    },
    "smoking-frozen-controls" = {
      frozen <- load_smoking_frozen_controls_for_hardening()
      list(
        data = frozen$data,
        args = list(
          y = "lcigsale",
          ivar = "state",
          tvar = "year",
          d = "d",
          post = "post",
          controls = frozen$controls
        )
      )
    },
    stop(sprintf("unknown hardening case: %s", case_id))
  )

  result <- withCallingHandlers(
    do.call(
      lwdid_sensitivity,
      c(
        list(
          data = case$data,
          type = "all",
          verbose = FALSE
        ),
        case$args
      )
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    y_sd = stats::sd(case$data[[case$args$y]], na.rm = TRUE),
    result = result,
    warnings = sort(unique(warnings_seen))
  )
}

test_that("E8-03 hardening audit oracle matches current comprehensive outputs", {
  skip_if_no_parity_fixtures()
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_hardening_audit_path(
    "20260324-qa-parity-e8-03-comprehensive-numeric-audit.json"
  )

  raw_csv_path <- resolve_smoking_raw_csv_for_hardening()
  fixture_csv <- resolve_hardening_audit_path("e8_03_smoking_frozen_controls_fixture.csv")
  skip_if(!file.exists(oracle_path), paste("missing oracle:", oracle_path))
  skip_if(!file.exists(raw_csv_path), paste("missing raw csv:", raw_csv_path))
  skip_if(!file.exists(fixture_csv), paste("missing fixture:", fixture_csv))

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)

  expect_identical(oracle$comparison$status, "passed")
  expect_identical(oracle$comparison$numeric_status, "passed")

  for (case in oracle$cases) {
    current <- run_hardening_case(case$id)
    current_audit <- compute_numeric_audit(
      current$result,
      y_sd = current$y_sd,
      att_bound_multiplier = case$att_bound_multiplier
    )

    expect_s3_class(current$result, "lwdid_sensitivity_comprehensive")
    .filter_glm_warnings <- function(w) {
      w[!grepl("^glm\\.fit:", w) & !grepl("^no non-missing arguments", w)]
    }
    expect_identical(
      sort(.filter_glm_warnings(current$warnings)),
      sort(.filter_glm_warnings(unlist(case$warnings, use.names = FALSE)))
    )
    expect_equal(current$y_sd, case$y_sd, tolerance = 1e-12)
    expect_equal(current_audit$numeric_field_count, as.integer(case$numeric_field_count))
    expect_identical(
      current_audit$nonfinite_paths,
      normalize_oracle_character(case$nonfinite_paths)
    )
    expect_identical(
      current_audit$out_of_range_paths,
      normalize_oracle_character(case$out_of_range_paths)
    )
    expect_equal(current_audit$max_att_ratio, case$max_att_ratio, tolerance = 1e-12)

    oracle_att_ratios <- normalize_oracle_att_ratios(case$att_ratio_by_path)
    expect_identical(sort(names(current_audit$att_ratio_by_path)), sort(names(oracle_att_ratios)))
    expect_equal(
      current_audit$att_ratio_by_path[names(oracle_att_ratios)],
      oracle_att_ratios,
      tolerance = 1e-12
    )
  }
})
