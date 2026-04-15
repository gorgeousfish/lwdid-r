local_suppress_lwdid_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (inherits(w, "lwdid_warning") && !inherits(w, "lwdid_data")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

make_wcb_fixture <- function() {
  utils::read.csv(
    "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e9_04_layer2_first_slice_fixture.csv",
    stringsAsFactors = FALSE
  )
}

make_wcb_smoking_small_g_fixture <- function() {
  utils::data("smoking", package = "lwdid", envir = environment())
  state_subset <- sort(unique(smoking$state))[1:10]
  smoking_small_g <- subset(smoking, state %in% state_subset)
  result <- suppressWarnings(
    lwdid(
      smoking_small_g,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      rolling = "demean",
      vce = "cluster",
      cluster_var = "state"
    )
  )

  data.frame(
    y_dot = result$.ri_data$y_trans,
    treated = result$.ri_data$d,
    state = state_subset
  )
}

make_wcb_ill_conditioned_fixture <- function() {
  d <- c(rep(1, 99), 1 + 1e-5)

  data.frame(
    Y = 3 + 2 * d + rep(1:10, each = 10),
    D = d,
    cluster = rep(1:10, each = 10)
  )
}

make_wcb_controls_fixture <- function() {
  data.frame(
    Y = c(1.30, 4.00, 1.95, 4.85, 0.85, 3.85, 1.65, 4.70, 1.05, 3.95),
    D = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    cluster = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    x1 = c(0.15, 0.85, 0.35, 1.05, 0.05, 0.95, 0.25, 1.10, 0.10, 0.90)
  )
}

make_wcb_g12_fixture <- function() {
  cluster <- rep(1:12, each = 3)
  treated <- rep(c(0, 0, 1), 12)
  cluster_level <- rep(seq(-1.1, 1.1, length.out = 12), each = 3)
  within_cluster <- rep(c(-0.25, 0.15, 0.35), 12)
  deterministic_noise <- c(
    rep(c(-0.10, 0.05, 0.20), 4),
    rep(c(0.08, -0.03, 0.12), 4),
    rep(c(-0.06, 0.09, 0.18), 4)
  )

  data.frame(
    Y = 2 + 0.6 * cluster_level + 1.4 * treated + within_cluster + deterministic_noise,
    D = treated,
    cluster = cluster
  )
}

make_wcb_story_local_shared_dgp <- function() {
  list(
    intercept = 10.0,
    cluster_effect_distribution = "normal(0, 2)",
    idiosyncratic_error_distribution = "normal(0, 1)",
    treatment_assignment = list(
      level = "cluster",
      default_probability = 0.5
    )
  )
}

make_wcb_story_local_decision_scenario <- function(
    scenario_id,
    true_tau,
    n_simulations,
    metric_name,
    decision_direction,
    obs_per_cluster = 20L,
    seed = 42L
) {
  list(
    id = scenario_id,
    seed = as.integer(seed),
    n_simulations = as.integer(n_simulations),
    true_tau = as.numeric(true_tau),
    G = 10L,
    obs_per_cluster = as.integer(obs_per_cluster),
    estimator = "wild_cluster_bootstrap",
    n_bootstrap = 199L,
    requested_n_bootstrap = 199L,
    weight_type = "rademacher",
    alpha = 0.05,
    impose_null = TRUE,
    metric_name = metric_name,
    decision_direction = decision_direction,
    pvalue_threshold = 0.05
  )
}

test_that("TC-9.4.1: .generate_weights returns a Rademacher draw", {
  set.seed(123)
  weights <- lwdid:::.generate_weights(8L, "rademacher")

  expect_type(weights, "double")
  expect_length(weights, 8L)
  expect_true(all(weights %in% c(-1, 1)))
})

test_that("TC-9.4.2: .generate_weights returns Mammen draws with the expected moments", {
  sqrt5 <- sqrt(5)
  phi <- (sqrt5 + 1) / 2
  support <- c(-(phi - 1), phi)

  set.seed(123)
  weights <- lwdid:::.generate_weights(20000L, "mammen")

  expect_type(weights, "double")
  expect_length(weights, 20000L)
  expect_true(all(weights %in% support))
  expect_equal(mean(weights), 0, tolerance = 0.03)
  expect_equal(mean(weights^2), 1, tolerance = 0.03)
  expect_equal(mean(weights^3), 1, tolerance = 0.08)
})

test_that("TC-9.4.3: .generate_weights returns Webb draws on the six-point support", {
  support <- c(-sqrt(3 / 2), -1, -sqrt(1 / 2), sqrt(1 / 2), 1, sqrt(3 / 2))

  set.seed(123)
  weights <- lwdid:::.generate_weights(60000L, "webb")
  weights_table <- table(factor(weights, levels = support))

  expect_type(weights, "double")
  expect_length(weights, 60000L)
  expect_equal(unname(sort(unique(weights))), unname(sort(support)), tolerance = 1e-12)
  expect_equal(as.numeric(weights_table) / length(weights), rep(1 / 6, 6), tolerance = 0.02)
  expect_equal(mean(weights), 0, tolerance = 0.02)
  expect_equal(mean(weights^2), 1, tolerance = 0.03)
})

test_that("TC-9.4.4: .generate_all_rademacher enumerates all 2^G sign patterns", {
  weights <- lwdid:::.generate_all_rademacher(5L)

  expect_equal(dim(weights), c(32L, 5L))
  expect_true(all(weights %in% c(-1, 1)))
  expect_equal(nrow(unique(weights)), 32L)
})

test_that("TC-9.4.5: impose_null changes the WCB boundary on the shared first-slice fixture", {
  fixture <- make_wcb_fixture()

  impose_null_true <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      use_fwildclusterboot = FALSE
    )
  )
  impose_null_false <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = FALSE,
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(impose_null_true, "lwdid_wcb_result")
  expect_s3_class(impose_null_false, "lwdid_wcb_result")
  expect_equal(impose_null_true$att, impose_null_false$att, tolerance = 1e-12)
  expect_equal(
    impose_null_true$original_se,
    impose_null_false$original_se,
    tolerance = 1e-12
  )
  expect_identical(impose_null_true$actual_n_bootstrap, 32L)
  expect_identical(impose_null_false$actual_n_bootstrap, 32L)
  expect_true(isTRUE(impose_null_true$full_enumeration))
  expect_true(isTRUE(impose_null_false$full_enumeration))
  expect_identical(impose_null_true$ci_method, "percentile_t")
  expect_identical(impose_null_false$ci_method, "percentile")
  expect_equal(impose_null_true$pvalue, 0.0625, tolerance = 1e-12)
  expect_equal(impose_null_false$pvalue, 0.9375, tolerance = 1e-12)
  expect_equal(impose_null_true$ci_lower, 1.0615000000000006, tolerance = 1e-12)
  expect_equal(impose_null_true$ci_upper, 3.7635000000000005, tolerance = 1e-12)
  expect_equal(impose_null_false$ci_lower, 2.217083333333334, tolerance = 1e-12)
  expect_equal(impose_null_false$ci_upper, 2.6079166666666667, tolerance = 1e-12)
  expect_gt(abs(impose_null_true$pvalue - impose_null_false$pvalue), 0.5)
  expect_gt(
    abs(
      (impose_null_true$ci_upper - impose_null_true$ci_lower) -
        (impose_null_false$ci_upper - impose_null_false$ci_lower)
    ),
    1
  )
})

test_that("TC-9.4.9: .fast_cluster_se matches the manual cluster sandwich formula", {
  fixture <- make_wcb_fixture()
  x <- cbind("(Intercept)" = 1, D = fixture$D)
  y <- fixture$Y
  xtx_inv <- solve(crossprod(x))
  beta_hat <- drop(xtx_inv %*% crossprod(x, y))
  residuals <- y - drop(x %*% beta_hat)
  cluster_index <- as.integer(as.factor(fixture$cluster))
  g <- length(unique(cluster_index))

  fast_se <- lwdid:::.fast_cluster_se(
    residuals = residuals,
    X = x,
    obs_cluster_idx = cluster_index,
    G = g,
    XtX_inv = xtx_inv,
    param_idx = 2L
  )

  cluster_scores <- lapply(split(seq_len(nrow(fixture)), cluster_index), function(idx) {
    drop(crossprod(x[idx, , drop = FALSE], residuals[idx]))
  })
  meat <- Reduce(`+`, lapply(cluster_scores, tcrossprod))
  correction <- (g / (g - 1)) * ((nrow(x) - 1) / (nrow(x) - ncol(x)))
  var_beta <- correction * xtx_inv %*% meat %*% xtx_inv
  expected_se <- sqrt(as.numeric(var_beta[2L, 2L]))

  expect_type(fast_se, "double")
  expect_equal(fast_se, expected_se, tolerance = 1e-12)
  expect_gt(fast_se, 0)
})

test_that("TC-9.4.12: ill-conditioned WCB design emits a numerical warning", {
  fixture <- make_wcb_ill_conditioned_fixture()
  design_cond <- kappa(crossprod(cbind(1, fixture$D)), exact = TRUE)

  expect_gt(design_cond, 1e10)

  result <- NULL
  expect_warning(
    result <- lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher"
    ),
    regexp = "condition number",
    class = "lwdid_numerical"
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_equal(result$n_clusters, 10L)
})

test_that("TC-9.4.11: restricted_model changes the controls-aware impose-null path", {
  fixture <- make_wcb_controls_fixture()

  intercept_only <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "intercept_only"
    )
  )

  with_controls <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "with_controls"
    )
  )

  expect_s3_class(intercept_only, "lwdid_wcb_result")
  expect_s3_class(with_controls, "lwdid_wcb_result")
  expect_equal(intercept_only$restricted_model, "intercept_only")
  expect_equal(with_controls$restricted_model, "with_controls")
  expect_equal(intercept_only$actual_n_bootstrap, 32L)
  expect_true(isTRUE(intercept_only$full_enumeration))
  expect_equal(with_controls$actual_n_bootstrap, 32L)
  expect_true(isTRUE(with_controls$full_enumeration))
  expect_equal(intercept_only$att, -0.05445544554455495, tolerance = 1e-12)
  expect_equal(intercept_only$original_se, with_controls$original_se, tolerance = 1e-12)
  expect_equal(intercept_only$pvalue, 0.8125, tolerance = 1e-12)
  expect_equal(with_controls$pvalue, 0.84375, tolerance = 1e-12)
  expect_equal(with_controls$att, intercept_only$att, tolerance = 1e-12)
  expect_gt(abs(with_controls$pvalue - intercept_only$pvalue), 1e-6)
  expect_gt(abs(with_controls$ci_lower - intercept_only$ci_lower), 1e-6)
  expect_gt(abs(with_controls$ci_upper - intercept_only$ci_upper), 1e-6)
})

test_that("controls-aware intercept_only bootstrap SE matches the Python exact-enumeration oracle", {
  fixture <- make_wcb_controls_fixture()

  intercept_only <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "intercept_only"
    )
  )

  expect_s3_class(intercept_only, "lwdid_wcb_result")
  expect_equal(intercept_only$se_bootstrap, 1.9685369632567724, tolerance = 1e-12)
  expect_equal(intercept_only$ci_lower, -0.8219716790616101, tolerance = 1e-12)
  expect_equal(intercept_only$ci_upper, 0.7130607879725002, tolerance = 1e-12)
  expect_equal(intercept_only$pvalue, 0.8125, tolerance = 1e-12)
})

test_that("native first-slice runner materializes the public WCB contract", {
  fixture <- make_wcb_fixture()

  result <- local_suppress_lwdid_warnings(
    lwdid:::.wcb_native(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      alpha = 0.05,
      seed = 1L,
      impose_null = TRUE,
      full_enumeration = NULL,
      restricted_model = "with_controls"
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_equal(result$method, "native")
  expect_equal(result$requested_n_bootstrap, 999L)
  expect_equal(result$actual_n_bootstrap, 32L)
  expect_true(isTRUE(result$full_enumeration))
  expect_equal(result$restricted_model, "with_controls")
})

test_that("no-controls native path materializes bootstrap dispersion outputs", {
  fixture <- make_wcb_fixture()

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 1L,
      impose_null = TRUE
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_equal(result$actual_n_bootstrap, 32L)
  expect_true(isTRUE(result$full_enumeration))
  expect_type(result$se_bootstrap, "double")
  expect_true(is.finite(result$se_bootstrap))
  expect_gt(result$se_bootstrap, 0)
  expect_length(result$t_stats_bootstrap, 32L)
  expect_true(all(is.finite(result$t_stats_bootstrap)))
})

test_that("native WCB result materializes bootstrap ATT support for diagnostics", {
  fixture <- make_wcb_controls_fixture()

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "with_controls"
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_type(result$att_bootstrap, "double")
  expect_length(result$att_bootstrap, 32L)
  expect_true(all(is.finite(result$att_bootstrap)))
  if (is.double(result$att_bootstrap) && length(result$att_bootstrap) == 32L) {
    expect_equal(
      result$se_bootstrap,
      sqrt(mean((result$att_bootstrap - mean(result$att_bootstrap))^2)),
      tolerance = 1e-12
    )
  }
})

test_that("native WCB result materializes CI critical values for diagnostics", {
  fixture <- make_wcb_controls_fixture()

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "with_controls"
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_type(result$ci_lower_critical, "double")
  expect_type(result$ci_upper_critical, "double")
  expect_length(result$ci_lower_critical, 1L)
  expect_length(result$ci_upper_critical, 1L)
  expect_true(is.finite(result$ci_lower_critical))
  expect_true(is.finite(result$ci_upper_critical))
  expect_equal(
    result$ci_lower_critical,
    (result$att - result$ci_lower) / result$se_bootstrap,
    tolerance = 1e-12
  )
  expect_equal(
    result$ci_upper_critical,
    (result$ci_upper - result$att) / result$se_bootstrap,
    tolerance = 1e-12
  )
})

test_that("native WCB result exposes scale-aware critical diagnostics", {
  fixture <- make_wcb_controls_fixture()

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "with_controls"
    )
  )

  valid_t <- abs(result$t_stats_bootstrap[is.finite(result$t_stats_bootstrap)])

  expect_s3_class(result, "lwdid_wcb_result")
  expect_identical(result$t_stat_critical_scale, "original_t")
  expect_identical(result$ci_critical_scale, "bootstrap_se")
  expect_equal(
    result$t_stat_critical,
    stats::quantile(valid_t, probs = 0.95, names = FALSE, type = 7),
    tolerance = 1e-12
  )
  expect_equal(
    result$bootstrap_scale_t_stat,
    abs(result$att / result$se_bootstrap),
    tolerance = 1e-12
  )
  expect_equal(
    result$rejection_rate,
    mean(valid_t >= abs(result$t_stat_original)),
    tolerance = 1e-12
  )
})

test_that("native WCB records percentile CI metadata when impose_null is FALSE", {
  fixture <- make_wcb_fixture()

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 1L,
      impose_null = FALSE
    )
  )

  att_valid <- result$att_bootstrap[is.finite(result$att_bootstrap)]

  expect_s3_class(result, "lwdid_wcb_result")
  expect_false(result$impose_null)
  expect_equal(result$ci_method, "percentile")
  expect_equal(
    result$ci_lower,
    as.numeric(stats::quantile(att_valid, probs = 0.025, names = FALSE, type = 7)),
    tolerance = 1e-12
  )
  expect_equal(
    result$ci_upper,
    as.numeric(stats::quantile(att_valid, probs = 0.975, names = FALSE, type = 7)),
    tolerance = 1e-12
  )
})

test_that("native WCB accepts non-rademacher weights without full enumeration", {
  fixture <- make_wcb_fixture()

  for (weight_type in c("mammen", "webb")) {
    result <- local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        n_bootstrap = 199L,
        weight_type = weight_type,
        seed = 42L,
        impose_null = TRUE
      )
    )

    expect_s3_class(result, "lwdid_wcb_result")
    expect_equal(result$weight_type, weight_type)
    expect_equal(result$requested_n_bootstrap, 199L)
    expect_equal(result$actual_n_bootstrap, 199L)
    expect_false(isTRUE(result$full_enumeration))
    expect_true(is.finite(result$att))
    expect_true(is.finite(result$original_se))
    expect_true(is.finite(result$pvalue))
  }
})

test_that("native WCB result materializes restricted null-base diagnostics", {
  fixture <- make_wcb_controls_fixture()
  control_matrix <- as.matrix(fixture["x1"])
  restricted_matrix <- cbind("(Intercept)" = 1, control_matrix)
  expected_beta <- drop(
    solve(crossprod(restricted_matrix), crossprod(restricted_matrix, fixture$Y))
  )
  expected_fitted <- drop(restricted_matrix %*% expected_beta)
  expected_resid <- fixture$Y - expected_fitted

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "with_controls"
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_type(result$restricted_coefficients, "double")
  expect_equal(unname(result$restricted_coefficients), unname(expected_beta), tolerance = 1e-12)
  expect_equal(unname(result$fitted_base), unname(expected_fitted), tolerance = 1e-12)
  expect_equal(unname(result$resid_base), unname(expected_resid), tolerance = 1e-12)
  expect_equal(
    unname(result$fitted_base + result$resid_base),
    unname(fixture$Y),
    tolerance = 1e-12
  )
})

test_that("native helper materializes exact-enumeration pseudo-outcome support", {
  fixture <- make_wcb_controls_fixture()
  cluster_index <- as.integer(as.factor(fixture$cluster))
  weight_info <- lwdid:::.resolve_wild_bootstrap_weights(
    n_clusters = length(unique(cluster_index)),
    requested_n_bootstrap = 999L,
    weight_type = "rademacher",
    full_enumeration = NULL
  )
  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "with_controls"
    )
  )

  expected_weights <- weight_info$weights[, cluster_index, drop = FALSE]
  expected <- matrix(
    result$fitted_base,
    nrow = nrow(expected_weights),
    ncol = ncol(expected_weights),
    byrow = TRUE
  ) + expected_weights * matrix(
    result$resid_base,
    nrow = nrow(expected_weights),
    ncol = ncol(expected_weights),
    byrow = TRUE
  )
  pseudo_outcomes <- lwdid:::.build_wcb_pseudo_outcomes(
    fitted_base = result$fitted_base,
    resid_base = result$resid_base,
    weights = weight_info$weights,
    cluster_index = cluster_index
  )

  expect_type(pseudo_outcomes, "double")
  expect_equal(dim(pseudo_outcomes), c(32L, nrow(fixture)))
  expect_equal(pseudo_outcomes, expected, tolerance = 1e-12)
})

test_that("TC-9.4.6 TC-9.4.7 TC-9.4.14 TC-9.4.17: public first slice is deterministic and exposes the minimal result contract", {
  fixture <- make_wcb_fixture()

  set.seed(1)
  result_seed_1 <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      controls = NULL
    )
  )

  set.seed(999)
  result_seed_999 <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      controls = NULL
    )
  )

  expect_s3_class(result_seed_1, "lwdid_wcb_result")
  expect_equal(result_seed_1$requested_n_bootstrap, 999L)
  expect_equal(result_seed_1$actual_n_bootstrap, 32L)
  expect_true(isTRUE(result_seed_1$full_enumeration))
  expect_equal(result_seed_1$weight_type, "rademacher")
  expect_equal(result_seed_1$rejection_rate, result_seed_1$pvalue)
  expect_gt(result_seed_1$pvalue, 0)

  expect_equal(result_seed_1$att, result_seed_999$att, tolerance = 1e-12)
  expect_equal(result_seed_1$ci_lower, result_seed_999$ci_lower, tolerance = 1e-12)
  expect_equal(result_seed_1$ci_upper, result_seed_999$ci_upper, tolerance = 1e-12)
  expect_equal(result_seed_1$pvalue, result_seed_999$pvalue, tolerance = 1e-12)
})

test_that("TC-9.4.20: degenerate single-cluster input returns a NaN result object", {
  degenerate <- data.frame(
    Y = rep(1, 6),
    D = c(0, 0, 0, 1, 1, 1),
    cluster = rep(1, 6)
  )

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = degenerate,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher"
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_true(is.nan(result$pvalue))
  expect_true(is.nan(result$rejection_rate))
  expect_equal(result$n_clusters, 1L)
})

test_that("TC-9.4.8 TC-9.4.22: fully degenerate input returns a NaN result that print and summary can render", {
  degenerate <- data.frame(
    Y = rep(1, 10),
    D = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    cluster = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  )

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = degenerate,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      seed = 1L,
      impose_null = TRUE
    )
  )

  print_output <- capture.output(print(result))
  summary_output <- capture.output(summary(result))

  expect_s3_class(result, "lwdid_wcb_result")
  expect_true(is.nan(result$original_se))
  expect_true(is.nan(result$pvalue))
  expect_true(is.nan(result$ci_lower))
  expect_true(is.nan(result$ci_upper))
  expect_true(any(grepl("Wild Cluster Bootstrap Result", print_output, fixed = TRUE)))
  expect_true(any(grepl("ATT:", print_output, fixed = TRUE)))
  expect_true(any(grepl("Wild Cluster Bootstrap Result", summary_output, fixed = TRUE)))
  expect_true(any(grepl("P-value:", summary_output, fixed = TRUE)))
})

test_that("Layer 3 smoking small-G replay stays deterministic on the public first slice", {
  fixture <- make_wcb_smoking_small_g_fixture()

  result_seed_123 <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "y_dot",
      d = "treated",
      cluster_var = "state",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      controls = NULL,
      seed = 123L
    )
  )
  result_seed_999 <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "y_dot",
      d = "treated",
      cluster_var = "state",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      controls = NULL,
      seed = 999L
    )
  )

  expect_s3_class(result_seed_123, "lwdid_wcb_result")
  expect_equal(result_seed_123$n_clusters, 10L)
  expect_equal(result_seed_123$requested_n_bootstrap, 999L)
  expect_equal(result_seed_123$actual_n_bootstrap, 1024L)
  expect_true(isTRUE(result_seed_123$full_enumeration))
  expect_equal(result_seed_123$weight_type, "rademacher")
  expect_equal(result_seed_123$att, -0.442085549845333, tolerance = 1e-12)
  expect_equal(result_seed_123$original_se, 0.0475307783927852, tolerance = 1e-12)
  expect_equal(result_seed_123$ci_lower, -0.979784651460434, tolerance = 1e-12)
  expect_equal(result_seed_123$ci_upper, 0.0956135517697672, tolerance = 1e-12)
  expect_equal(result_seed_123$pvalue, 0.1796875, tolerance = 1e-12)
  expect_equal(result_seed_123$rejection_rate, result_seed_123$pvalue, tolerance = 1e-12)

  expect_equal(result_seed_123$att, result_seed_999$att, tolerance = 1e-12)
  expect_equal(result_seed_123$original_se, result_seed_999$original_se, tolerance = 1e-12)
  expect_equal(result_seed_123$ci_lower, result_seed_999$ci_lower, tolerance = 1e-12)
  expect_equal(result_seed_123$ci_upper, result_seed_999$ci_upper, tolerance = 1e-12)
  expect_equal(result_seed_123$pvalue, result_seed_999$pvalue, tolerance = 1e-12)
})

test_that("TC-9.4.10: hard-coded Python parity snapshot matches the native exact-enumeration contract", {
  fixture <- make_wcb_fixture()

  expected <- list(
    att = 2.4125000000000005,
    original_se = 0.115463789016618,
    ci_lower = 1.0615000000000006,
    ci_upper = 3.7635000000000005,
    pvalue = 0.0625,
    actual_n_bootstrap = 32L,
    weight_type = "rademacher",
    ci_method = "percentile_t"
  )

  run_contract <- function(seed, requested_n_bootstrap) {
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        n_bootstrap = requested_n_bootstrap,
        weight_type = "rademacher",
        seed = seed,
        impose_null = TRUE,
        use_fwildclusterboot = FALSE
      )
    )
  }

  seed_123_requested_199 <- run_contract(seed = 123L, requested_n_bootstrap = 199L)
  seed_999_requested_199 <- run_contract(seed = 999L, requested_n_bootstrap = 199L)
  seed_123_requested_999 <- run_contract(seed = 123L, requested_n_bootstrap = 999L)

  for (result in list(
    seed_123_requested_199,
    seed_999_requested_199,
    seed_123_requested_999
  )) {
    expect_s3_class(result, "lwdid_wcb_result")
    expect_equal(result$att, expected$att, tolerance = 1e-12)
    expect_equal(result$original_se, expected$original_se, tolerance = 1e-12)
    expect_equal(result$ci_lower, expected$ci_lower, tolerance = 1e-12)
    expect_equal(result$ci_upper, expected$ci_upper, tolerance = 1e-12)
    expect_equal(result$pvalue, expected$pvalue, tolerance = 1e-12)
    expect_equal(result$actual_n_bootstrap, expected$actual_n_bootstrap)
    expect_true(isTRUE(result$full_enumeration))
    expect_equal(result$weight_type, expected$weight_type)
    expect_equal(result$ci_method, expected$ci_method)
  }

  expect_equal(seed_123_requested_199$requested_n_bootstrap, 199L)
  expect_equal(seed_999_requested_199$requested_n_bootstrap, 199L)
  expect_equal(seed_123_requested_999$requested_n_bootstrap, 999L)

  expect_equal(seed_123_requested_199$att, seed_999_requested_199$att, tolerance = 1e-12)
  expect_equal(seed_123_requested_199$original_se, seed_999_requested_199$original_se, tolerance = 1e-12)
  expect_equal(seed_123_requested_199$ci_lower, seed_999_requested_199$ci_lower, tolerance = 1e-12)
  expect_equal(seed_123_requested_199$ci_upper, seed_999_requested_199$ci_upper, tolerance = 1e-12)
  expect_equal(seed_123_requested_199$pvalue, seed_999_requested_199$pvalue, tolerance = 1e-12)

  expect_equal(seed_123_requested_199$att, seed_123_requested_999$att, tolerance = 1e-12)
  expect_equal(seed_123_requested_199$original_se, seed_123_requested_999$original_se, tolerance = 1e-12)
  expect_equal(seed_123_requested_199$ci_lower, seed_123_requested_999$ci_lower, tolerance = 1e-12)
  expect_equal(seed_123_requested_199$ci_upper, seed_123_requested_999$ci_upper, tolerance = 1e-12)
  expect_equal(seed_123_requested_199$pvalue, seed_123_requested_999$pvalue, tolerance = 1e-12)
})

test_that("TC-9.4.13: .wcb_via_fwildclusterboot maps the boottest contract into lwdid_wcb_result", {
  fixture <- make_wcb_controls_fixture()
  recorded <- new.env(parent = emptyenv())

  fake_boottest <- function(
      object, param, clustid, B, type, impose_null, seed,
      conf_int = TRUE, sign_level = 0.05, ...
  ) {
    recorded$formula <- stats::formula(object)
    recorded$param <- param
    recorded$clustid <- clustid
    recorded$B <- B
    recorded$type <- type
    recorded$impose_null <- impose_null
    recorded$seed <- seed
    recorded$conf_int <- conf_int
    recorded$sign_level <- sign_level

    list(
      p_val = 0.125,
      conf_int = c(-0.5, 1.5),
      t_stat = 2.75,
      boot_iter = B,
      t_boot = numeric(0)
    )
  }

  result <- lwdid:::.wcb_via_fwildclusterboot(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    n_bootstrap = 399L,
    weight_type = "webb",
    alpha = 0.10,
    seed = 99L,
    impose_null = FALSE,
    boottest_fn = fake_boottest
  )

  expect_equal(deparse(recorded$formula), "Y ~ D + x1")
  expect_equal(recorded$param, "D")
  expect_equal(recorded$clustid, "cluster")
  expect_equal(recorded$B, 399L)
  expect_equal(recorded$type, "webb")
  expect_false(recorded$impose_null)
  expect_equal(recorded$seed, 99L)
  expect_true(recorded$conf_int)
  expect_equal(recorded$sign_level, 0.10, tolerance = 1e-12)

  expect_s3_class(result, "lwdid_wcb_result")
  expect_equal(result$method, "fwildclusterboot")
  expect_true(isTRUE(result$use_fwildclusterboot))
  expect_equal(result$weight_type, "webb")
  expect_equal(result$requested_n_bootstrap, 399L)
  expect_equal(result$actual_n_bootstrap, 399L)
  expect_false(isTRUE(result$full_enumeration))
  expect_false(result$impose_null)
  expect_equal(result$ci_method, "percentile_t")
  expect_equal(result$pvalue, 0.125, tolerance = 1e-12)
  expect_equal(result$ci_lower, -0.5, tolerance = 1e-12)
  expect_equal(result$ci_upper, 1.5, tolerance = 1e-12)
  expect_equal(result$t_stat_original, 2.75, tolerance = 1e-12)
  expect_equal(result$att, unname(stats::coef(stats::lm(Y ~ D + x1, data = fixture))[["D"]]), tolerance = 1e-12)
})

test_that("TC-9.4.13: .wcb_via_fwildclusterboot retries boottest interfaces that do not accept seed", {
  fixture <- make_wcb_controls_fixture()
  recorded <- new.env(parent = emptyenv())
  recorded$draws <- numeric(0)

  fake_boottest <- function(
      object, param, clustid, B, type, impose_null,
      conf_int = TRUE, sign_level = 0.05
  ) {
    recorded$formula <- stats::formula(object)
    recorded$param <- param
    recorded$clustid <- clustid
    recorded$B <- B
    recorded$type <- type
    recorded$impose_null <- impose_null
    recorded$conf_int <- conf_int
    recorded$sign_level <- sign_level
    recorded$draws <- c(recorded$draws, stats::runif(1))

    list(
      p_val = 0.2,
      conf_int = c(-0.25, 0.75),
      t_stat = 1.5,
      boot_iter = B,
      t_boot = numeric(0)
    )
  }

  result_seed_11_a <- lwdid:::.wcb_via_fwildclusterboot(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    n_bootstrap = 399L,
    weight_type = "webb",
    alpha = 0.10,
    seed = 11L,
    impose_null = FALSE,
    boottest_fn = fake_boottest
  )
  result_seed_11_b <- lwdid:::.wcb_via_fwildclusterboot(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    n_bootstrap = 399L,
    weight_type = "webb",
    alpha = 0.10,
    seed = 11L,
    impose_null = FALSE,
    boottest_fn = fake_boottest
  )
  result_seed_12 <- lwdid:::.wcb_via_fwildclusterboot(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    n_bootstrap = 399L,
    weight_type = "webb",
    alpha = 0.10,
    seed = 12L,
    impose_null = FALSE,
    boottest_fn = fake_boottest
  )

  expect_equal(deparse(recorded$formula), "Y ~ D + x1")
  expect_equal(recorded$param, "D")
  expect_equal(recorded$clustid, "cluster")
  expect_equal(recorded$B, 399L)
  expect_equal(recorded$type, "webb")
  expect_false(recorded$impose_null)
  expect_true(recorded$conf_int)
  expect_equal(recorded$sign_level, 0.10, tolerance = 1e-12)

  expect_length(recorded$draws, 3L)
  expect_equal(recorded$draws[[1L]], recorded$draws[[2L]], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(recorded$draws[[1L]], recorded$draws[[3L]])))

  expect_s3_class(result_seed_11_a, "lwdid_wcb_result")
  expect_equal(result_seed_11_a$method, "fwildclusterboot")
  expect_true(isTRUE(result_seed_11_a$use_fwildclusterboot))
  expect_equal(result_seed_11_a$pvalue, 0.2, tolerance = 1e-12)
  expect_equal(result_seed_11_a$ci_lower, -0.25, tolerance = 1e-12)
  expect_equal(result_seed_11_a$ci_upper, 0.75, tolerance = 1e-12)
  expect_equal(result_seed_11_b$pvalue, result_seed_11_a$pvalue, tolerance = 1e-12)
  expect_equal(result_seed_12$pvalue, result_seed_11_a$pvalue, tolerance = 1e-12)
})

test_that("TC-9.4.13: .wcb_via_fwildclusterboot keeps backend B above the fwildclusterboot minimum under full enumeration", {
  fixture <- make_wcb_fixture()
  recorded <- new.env(parent = emptyenv())

  fake_boottest <- function(
      object, param, clustid, B, type, impose_null, seed,
      conf_int = TRUE, sign_level = 0.05, ...
  ) {
    recorded$B <- B
    recorded$type <- type
    recorded$impose_null <- impose_null
    recorded$seed <- seed

    list(
      p_val = 0.125,
      conf_int = c(-0.5, 1.5),
      t_stat = 2.75,
      boot_iter = B,
      t_boot = numeric(0)
    )
  }

  result <- lwdid:::.wcb_via_fwildclusterboot(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = NULL,
    n_bootstrap = 199L,
    requested_n_bootstrap = 199L,
    weight_type = "rademacher",
    alpha = 0.05,
    seed = 1L,
    impose_null = TRUE,
    full_enumeration = TRUE,
    boottest_fn = fake_boottest
  )

  expect_identical(recorded$B, 199L)
  expect_identical(recorded$type, "rademacher")
  expect_true(recorded$impose_null)
  expect_identical(recorded$seed, 1L)
  expect_true(isTRUE(result$full_enumeration))
  expect_identical(result$actual_n_bootstrap, 32L)
})

test_that("TC-9.4.13: wild_cluster_bootstrap dispatches to fwildclusterboot when the adapter is available", {
  fixture <- make_wcb_controls_fixture()
  sentinel <- structure(
    list(
      att = 0.25,
      original_se = 0.1,
      restricted_coefficients = NULL,
      fitted_base = numeric(0),
      resid_base = numeric(0),
      se_bootstrap = NaN,
      ci_lower = -0.5,
      ci_upper = 1.5,
      ci_lower_critical = NaN,
      ci_upper_critical = NaN,
      pvalue = 0.125,
      rejection_rate = 0.125,
      n_clusters = 5L,
      requested_n_bootstrap = 399L,
      actual_n_bootstrap = 399L,
      n_bootstrap = 399L,
      weight_type = "webb",
      t_stat_original = 2.75,
      att_bootstrap = numeric(0),
      t_stats_bootstrap = numeric(0),
      method = "fwildclusterboot",
      impose_null = FALSE,
      full_enumeration = FALSE,
      ci_method = "percentile_t",
      restricted_model = "with_controls",
      use_fwildclusterboot = TRUE
    ),
    class = "lwdid_wcb_result"
  )

  testthat::local_mocked_bindings(
    .fwildclusterboot_available = function() TRUE,
    .wcb_via_fwildclusterboot = function(...) sentinel,
    .package = "lwdid"
  )

  expect_no_warning(
    result <- lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 399L,
      weight_type = "webb",
      alpha = 0.10,
      seed = 99L,
      impose_null = FALSE,
      use_fwildclusterboot = TRUE
    )
  )

  expect_identical(result, sentinel)
})

test_that("TC-9.4.13: fwildclusterboot dispatch preserves rademacher full-enumeration semantics", {
  fixture <- make_wcb_fixture()
  recorded <- new.env(parent = emptyenv())
  sentinel <- structure(
    list(
      att = 0.25,
      original_se = 0.1,
      restricted_coefficients = NULL,
      fitted_base = numeric(0),
      resid_base = numeric(0),
      se_bootstrap = NaN,
      ci_lower = -0.5,
      ci_upper = 1.5,
      ci_lower_critical = NaN,
      ci_upper_critical = NaN,
      pvalue = 0.125,
      rejection_rate = 0.125,
      n_clusters = 5L,
      requested_n_bootstrap = 199L,
      actual_n_bootstrap = 32L,
      n_bootstrap = 32L,
      weight_type = "rademacher",
      t_stat_original = 2.75,
      att_bootstrap = numeric(0),
      t_stats_bootstrap = numeric(0),
      method = "fwildclusterboot",
      impose_null = TRUE,
      full_enumeration = TRUE,
      ci_method = "percentile_t",
      restricted_model = "with_controls",
      use_fwildclusterboot = TRUE
    ),
    class = "lwdid_wcb_result"
  )

  testthat::local_mocked_bindings(
    .fwildclusterboot_available = function() TRUE,
    .wcb_via_fwildclusterboot = function(...) {
      recorded$args <- list(...)
      sentinel
    },
    .package = "lwdid"
  )

  expect_no_warning(
    result <- lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 199L,
      weight_type = "rademacher",
      seed = 99L,
      use_fwildclusterboot = TRUE
    )
  )

  expect_identical(result, sentinel)
  expect_identical(recorded$args$requested_n_bootstrap, 199L)
  expect_identical(recorded$args$actual_n_bootstrap, 32L)
  expect_true(isTRUE(recorded$args$full_enumeration))
})

test_that("TC-9.4.13: .wcb_via_fwildclusterboot lifts tiny full-enumeration B to 100", {
  fixture <- make_wcb_fixture()
  recorded <- new.env(parent = emptyenv())

  fake_boottest <- function(
      object, param, clustid, B, type, impose_null, seed,
      conf_int = TRUE, sign_level = 0.05, ...
  ) {
    recorded$B <- B
    if (B < 100L) {
      stop(
        paste(
          "The function argument B is smaller than 100.",
          "The number of bootstrap iterations needs to be 100 or higher in",
          "order to guarantee that the root finding procedure used to find",
          "the confidence set works properly."
        ),
        call. = FALSE
      )
    }

    list(
      p_val = 0.125,
      conf_int = c(-0.5, 1.5),
      t_stat = 2.75,
      boot_iter = B,
      t_boot = numeric(0)
    )
  }

  result <- lwdid:::.wcb_via_fwildclusterboot(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    n_bootstrap = 31L,
    requested_n_bootstrap = 31L,
    weight_type = "rademacher",
    alpha = 0.05,
    seed = 1L,
    impose_null = TRUE,
    full_enumeration = TRUE,
    boottest_fn = fake_boottest
  )

  expect_identical(recorded$B, 100L)
  expect_identical(result$requested_n_bootstrap, 31L)
  expect_identical(result$actual_n_bootstrap, 32L)
  expect_true(isTRUE(result$full_enumeration))
  expect_identical(result$method, "fwildclusterboot")
})

test_that("TC-9.4.13: .wcb_via_fwildclusterboot refreshes B when actual_n_bootstrap is overridden", {
  fixture <- make_wcb_fixture()
  recorded <- new.env(parent = emptyenv())

  fake_boottest <- function(
      object, param, clustid, B, type, impose_null, seed,
      conf_int = TRUE, sign_level = 0.05, ...
  ) {
    recorded$B <- B
    list(
      p_val = 0.125,
      conf_int = c(-0.5, 1.5),
      t_stat = 2.75,
      boot_iter = B,
      t_boot = numeric(0)
    )
  }

  result <- lwdid:::.wcb_via_fwildclusterboot(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    n_bootstrap = 199L,
    requested_n_bootstrap = 199L,
    actual_n_bootstrap = 1024L,
    weight_type = "rademacher",
    alpha = 0.05,
    seed = 1L,
    impose_null = TRUE,
    full_enumeration = TRUE,
    boottest_fn = fake_boottest
  )

  expect_identical(recorded$B, 1024L)
  expect_identical(result$requested_n_bootstrap, 199L)
  expect_identical(result$actual_n_bootstrap, 1024L)
  expect_true(isTRUE(result$full_enumeration))
  expect_identical(result$method, "fwildclusterboot")
})

test_that("TC-9.4.13: .wcb_via_fwildclusterboot keeps successful backend calls package-silent", {
  fixture <- make_wcb_controls_fixture()

  noisy_boottest <- function(
      object, param, clustid, B, type, impose_null, seed,
      conf_int = TRUE, sign_level = 0.05, ...
  ) {
    message("simulated fwildclusterboot startup chatter")
    warning("simulated fwildclusterboot informational warning", call. = FALSE)

    list(
      p_val = 0.02,
      conf_int = c(1.8, 2.1),
      t_stat = 4.2,
      t_boot = c(4.1, 4.3)
    )
  }

  expect_no_message(
    expect_no_warning(
      result <- lwdid:::.wcb_via_fwildclusterboot(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 199L,
        weight_type = "rademacher",
        alpha = 0.05,
        seed = 1L,
        impose_null = TRUE,
        full_enumeration = TRUE,
        boottest_fn = noisy_boottest
      )
    )
  )

  expect_identical(result$method, "fwildclusterboot")
  expect_true(isTRUE(result$use_fwildclusterboot))
  expect_equal(result$ci_lower, 1.8, tolerance = 1e-12)
  expect_equal(result$ci_upper, 2.1, tolerance = 1e-12)
  expect_equal(result$pvalue, 0.02, tolerance = 1e-12)
})

test_that("TC-9.4.13: fwildclusterboot runtime failures warn and fall back to the native path", {
  fixture <- make_wcb_controls_fixture()

  testthat::local_mocked_bindings(
    .fwildclusterboot_available = function() TRUE,
    .wcb_via_fwildclusterboot = function(...) {
      stop("simulated fwildclusterboot adapter failure", call. = FALSE)
    },
    .package = "lwdid"
  )

  native_result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 399L,
      weight_type = "webb",
      alpha = 0.10,
      seed = 99L,
      impose_null = FALSE,
      use_fwildclusterboot = FALSE
    )
  )

  expect_warning(
    fallback_result <- local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    regexp = "fwildclusterboot",
    class = "lwdid_data"
  )

  expect_s3_class(fallback_result, "lwdid_wcb_result")
  expect_identical(fallback_result$method, "native")
  expect_false(isTRUE(fallback_result$use_fwildclusterboot))
  expect_equal(fallback_result$att, native_result$att, tolerance = 1e-12)
  expect_equal(fallback_result$original_se, native_result$original_se, tolerance = 1e-12)
  expect_equal(fallback_result$se_bootstrap, native_result$se_bootstrap, tolerance = 1e-12)
  expect_equal(fallback_result$pvalue, native_result$pvalue, tolerance = 1e-12)
  expect_equal(fallback_result$ci_lower, native_result$ci_lower, tolerance = 1e-12)
  expect_equal(fallback_result$ci_upper, native_result$ci_upper, tolerance = 1e-12)
  expect_identical(
    fallback_result$requested_n_bootstrap,
    native_result$requested_n_bootstrap
  )
  expect_identical(
    fallback_result$actual_n_bootstrap,
    native_result$actual_n_bootstrap
  )
})

test_that("TC-9.4.26: fwildclusterboot-unavailable warning carries fallback metadata", {
  fixture <- make_wcb_fixture()
  captured_warning <- NULL
  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  testthat::local_mocked_bindings(
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        n_bootstrap = 199L,
        weight_type = "rademacher",
        seed = 1L,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(captured_warning$detail, "fwildclusterboot_unavailable")
  expect_identical(
    captured_warning$action_taken,
    "falling back to native bootstrap path"
  )
  expect_identical(captured_warning$requested_n_bootstrap, 199L)
  expect_identical(captured_warning$actual_n_bootstrap, 32L)
  expect_true(isTRUE(captured_warning$full_enumeration))
  expect_identical(captured_warning$weight_type, "rademacher")
  expect_false(isTRUE(captured_warning$adapter_available))
  expect_null(captured_warning$source_error)
  expect_identical(
    captured_warning$blocker_boundary,
    runtime_diagnostics$blocker_boundary
  )
  expect_identical(
    captured_warning$toolchain_mismatch_detected,
    runtime_diagnostics$toolchain_mismatch_detected
  )
  expect_identical(
    captured_warning$makeconf_path,
    runtime_diagnostics$makeconf_path
  )
  expect_identical(
    captured_warning$makeconf_exists,
    runtime_diagnostics$makeconf_exists
  )
  expect_identical(
    captured_warning$makeconf_flibs,
    runtime_diagnostics$makeconf_flibs
  )
  expect_identical(
    captured_warning$makeconf_fc,
    runtime_diagnostics$makeconf_fc
  )
  expect_identical(
    captured_warning$makeconf_f77,
    runtime_diagnostics$makeconf_f77
  )
  expect_identical(
    captured_warning$makeconf_uses_stale_opt_gfortran,
    runtime_diagnostics$makeconf_uses_stale_opt_gfortran
  )
  expect_identical(
    captured_warning$runtime_hint,
    runtime_diagnostics$runtime_hint
  )
})

test_that("TC-9.4.26: runtime diagnostics expose Makeconf wiring for stale FLIBS paths", {
  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/homebrew/bin/gfortran",
        "F77 = /opt/homebrew/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_identical(
    runtime_diagnostics$makeconf_path,
    "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
  )
  expect_true(isTRUE(runtime_diagnostics$makeconf_exists))
  expect_match(runtime_diagnostics$makeconf_flibs, "/opt/gfortran", fixed = TRUE)
  expect_identical(runtime_diagnostics$makeconf_fc, "/opt/homebrew/bin/gfortran")
  expect_identical(runtime_diagnostics$makeconf_f77, "/opt/homebrew/bin/gfortran")
  expect_true(isTRUE(runtime_diagnostics$makeconf_uses_stale_opt_gfortran))
  expect_true(isTRUE(runtime_diagnostics$toolchain_mismatch_detected))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
})

test_that("TC-9.4.26: runtime diagnostics lift the blocker when user Makevars clears stale FLIBS", {
  makevars_user_path <- tempfile("wcb-makevars-", fileext = ".mk")
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    makevars_user_path,
    useBytes = TRUE
  )
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(makevars_user_path)
  }, add = TRUE)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_true(isTRUE(runtime_diagnostics$toolchain_mismatch_detected))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "r-universe-index-visible / r-universe-source-archive-403",
    fixed = TRUE
  )
})

test_that("TC-9.4.26: package-source-closure runtime diagnostics expose dependency status", {
  makevars_user_path <- tempfile("wcb-makevars-", fileext = ".mk")
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    makevars_user_path,
    useBytes = TRUE
  )
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(makevars_user_path)
  }, add = TRUE)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )
  expect_false(isTRUE(runtime_diagnostics$fwildclusterboot_available))
  expect_named(
    runtime_diagnostics$dependency_status,
    c("fwildclusterboot", "summclust", "JuliaConnectoR", "sitmo", "dqrng")
  )
  expect_true(all(vapply(
    runtime_diagnostics$dependency_status,
    function(value) is.logical(value) && length(value) == 1L,
    logical(1)
  )))
  expect_identical(
    runtime_diagnostics$dependency_status$fwildclusterboot,
    runtime_diagnostics$fwildclusterboot_available
  )
  expected_missing <- names(Filter(
    function(value) identical(value, FALSE),
    runtime_diagnostics$dependency_status
  ))
  expect_identical(
    runtime_diagnostics$missing_dependency_names,
    expected_missing
  )
  expect_identical(
    runtime_diagnostics$missing_dependency_count,
    as.integer(length(expected_missing))
  )
  expect_identical(
    runtime_diagnostics$direct_install_failure_node,
    "r-universe-index-visible-source-archive-403"
  )
  expect_identical(
    runtime_diagnostics$direct_install_probe_provider,
    "r-universe"
  )
  expect_true(isTRUE(
    runtime_diagnostics$direct_install_provider_index_visible
  ))
  expect_true(isTRUE(
    runtime_diagnostics$direct_install_source_archive_forbidden_detected
  ))
  expect_identical(
    runtime_diagnostics$direct_install_source_archive_status_code,
    403L
  )
  expect_false(isTRUE(runtime_diagnostics$github_direct_install_tested))
  expect_match(
    runtime_diagnostics$runtime_hint,
    "r-universe-index-visible / r-universe-source-archive-403",
    fixed = TRUE
  )
})

test_that("TC-9.4.26: Makevars override also lifts stale Makeconf wiring when current FLIBS no longer points at missing paths", {
  makevars_user_path <- tempfile("wcb-makevars-", fileext = ".mk")
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    makevars_user_path,
    useBytes = TRUE
  )
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(makevars_user_path)
  }, add = TRUE)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_identical(runtime_diagnostics$flibs_missing_paths, character(0))
  expect_true(isTRUE(runtime_diagnostics$makeconf_uses_stale_opt_gfortran))
  expect_true(isTRUE(runtime_diagnostics$toolchain_mismatch_detected))
  expect_true(isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "r-universe-index-visible / r-universe-source-archive-403",
    fixed = TRUE
  )
})

test_that("TC-9.4.26: runtime diagnostics discover the default ~/.R/Makevars when R_MAKEVARS_USER is unset", {
  home_dir <- tempfile("wcb-home-")
  makevars_dir <- file.path(home_dir, ".R")
  makevars_user_path <- file.path(makevars_dir, "Makevars")
  dir.create(makevars_dir, recursive = TRUE)
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    makevars_user_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.unsetenv("R_MAKEVARS_USER")

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_true(isTRUE(runtime_diagnostics$makevars_user_exists))
  expect_identical(
    runtime_diagnostics$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(runtime_diagnostics$makevars_user_source, "default-home")
  expect_false(isTRUE(runtime_diagnostics$makevars_user_file_empty))
  expect_true(isTRUE(runtime_diagnostics$makevars_user_has_toolchain_override))
  expect_true(isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )
})

test_that("TC-9.4.26: whitespace-only R_MAKEVARS_USER falls back to the default ~/.R/Makevars", {
  home_dir <- tempfile("wcb-home-blank-")
  makevars_dir <- file.path(home_dir, ".R")
  makevars_user_path <- file.path(makevars_dir, "Makevars")
  dir.create(makevars_dir, recursive = TRUE)
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    makevars_user_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.setenv(R_MAKEVARS_USER = "   ")

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_true(isTRUE(runtime_diagnostics$makevars_user_exists))
  expect_identical(
    runtime_diagnostics$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(runtime_diagnostics$makevars_user_source, "default-home")
  expect_false(isTRUE(runtime_diagnostics$makevars_user_file_empty))
  expect_true(isTRUE(runtime_diagnostics$makevars_user_has_toolchain_override))
  expect_true(isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )
})

test_that("TC-9.4.26: fwildclusterboot warnings ignore blank R_MAKEVARS_USER and surface default-home diagnostics", {
  home_dir <- tempfile("wcb-home-warning-")
  makevars_dir <- file.path(home_dir, ".R")
  makevars_user_path <- file.path(makevars_dir, "Makevars")
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  dir.create(makevars_dir, recursive = TRUE)
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    makevars_user_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.setenv(R_MAKEVARS_USER = "   ")

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(
    captured_warning$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_true(isTRUE(captured_warning$makevars_user_exists))
  expect_identical(captured_warning$makevars_user_source, "default-home")
  expect_false(isTRUE(captured_warning$makevars_user_file_empty))
  expect_true(isTRUE(captured_warning$makevars_user_has_toolchain_override))
  expect_true(isTRUE(captured_warning$makevars_override_clears_stale_flibs))
  expect_identical(
    captured_warning$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )
})

test_that("TC-9.4.26: empty default ~/.R/Makevars narrows the runtime blocker to a missing override", {
  home_dir <- tempfile("wcb-home-empty-")
  makevars_dir <- file.path(home_dir, ".R")
  makevars_user_path <- file.path(makevars_dir, "Makevars")
  dir.create(makevars_dir, recursive = TRUE)
  writeLines(character(0), makevars_user_path, useBytes = TRUE)

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.unsetenv("R_MAKEVARS_USER")

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_true(isTRUE(runtime_diagnostics$makevars_user_exists))
  expect_identical(
    runtime_diagnostics$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(runtime_diagnostics$makevars_user_source, "default-home")
  expect_true(isTRUE(runtime_diagnostics$makevars_user_file_empty))
  expect_false(isTRUE(runtime_diagnostics$makevars_user_has_toolchain_override))
  expect_false(isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "~/.R/Makevars file is present but empty",
    fixed = TRUE
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "R_MAKEVARS_USER",
    fixed = TRUE
  )
})

test_that("TC-9.4.26: non-empty default ~/.R/Makevars gets a dedicated no-override hint", {
  home_dir <- tempfile("wcb-home-nonempty-default-")
  makevars_dir <- file.path(home_dir, ".R")
  makevars_user_path <- file.path(makevars_dir, "Makevars")
  dir.create(makevars_dir, recursive = TRUE)
  writeLines(
    c(
      "CFLAGS = -O2",
      "PKG_CPPFLAGS = -DUSE_DEFAULT_HOME_PROBE"
    ),
    makevars_user_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.unsetenv("R_MAKEVARS_USER")

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_true(isTRUE(runtime_diagnostics$makevars_user_exists))
  expect_identical(
    runtime_diagnostics$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(runtime_diagnostics$makevars_user_source, "default-home")
  expect_false(isTRUE(runtime_diagnostics$makevars_user_file_empty))
  expect_false(isTRUE(runtime_diagnostics$makevars_user_explicit_missing))
  expect_false(isTRUE(runtime_diagnostics$makevars_user_has_toolchain_override))
  expect_false(isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "default ~/.R/Makevars file is non-empty",
    fixed = TRUE
  )
  expect_false(grepl(
    "discovered user Makevars file",
    runtime_diagnostics$runtime_hint,
    fixed = TRUE
  ))
})

test_that("TC-9.4.26: fwildclusterboot warnings surface non-empty default ~/.R/Makevars diagnostics", {
  home_dir <- tempfile("wcb-home-nonempty-default-warning-")
  makevars_dir <- file.path(home_dir, ".R")
  makevars_user_path <- file.path(makevars_dir, "Makevars")
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  dir.create(makevars_dir, recursive = TRUE)
  writeLines(
    c(
      "CFLAGS = -O2",
      "PKG_CPPFLAGS = -DUSE_DEFAULT_HOME_WARNING_PROBE"
    ),
    makevars_user_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.unsetenv("R_MAKEVARS_USER")

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(
    captured_warning$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_true(isTRUE(captured_warning$makevars_user_exists))
  expect_identical(captured_warning$makevars_user_source, "default-home")
  expect_false(isTRUE(captured_warning$makevars_user_file_empty))
  expect_false(isTRUE(captured_warning$makevars_user_explicit_missing))
  expect_false(isTRUE(captured_warning$makevars_user_has_toolchain_override))
  expect_false(isTRUE(captured_warning$makevars_override_clears_stale_flibs))
  expect_identical(
    captured_warning$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    captured_warning$runtime_hint,
    "default ~/.R/Makevars file is non-empty",
    fixed = TRUE
  )
  expect_false(grepl(
    "discovered user Makevars file",
    captured_warning$runtime_hint,
    fixed = TRUE
  ))
})

test_that("TC-9.4.26: empty explicit R_MAKEVARS_USER gets a dedicated empty-override hint", {
  home_dir <- tempfile("wcb-home-empty-explicit-override-")
  makevars_user_path <- file.path(home_dir, "override.mk")
  dir.create(home_dir, recursive = TRUE)
  writeLines(character(0), makevars_user_path, useBytes = TRUE)

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_identical(
    runtime_diagnostics$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(runtime_diagnostics$makevars_user_source, "env-var")
  expect_true(isTRUE(runtime_diagnostics$makevars_user_exists))
  expect_true(isTRUE(runtime_diagnostics$makevars_user_file_empty))
  expect_false(isTRUE(runtime_diagnostics$makevars_user_explicit_missing))
  expect_false(isTRUE(runtime_diagnostics$makevars_user_has_toolchain_override))
  expect_false(isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "configured R_MAKEVARS_USER file is present but empty",
    fixed = TRUE
  )
  expect_false(grepl(
    "discovered user Makevars file",
    runtime_diagnostics$runtime_hint,
    fixed = TRUE
  ))
})

test_that("TC-9.4.26: fwildclusterboot warnings surface empty explicit R_MAKEVARS_USER diagnostics", {
  home_dir <- tempfile("wcb-home-empty-explicit-override-warning-")
  makevars_user_path <- file.path(home_dir, "override.mk")
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  dir.create(home_dir, recursive = TRUE)
  writeLines(character(0), makevars_user_path, useBytes = TRUE)

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(
    captured_warning$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(captured_warning$makevars_user_source, "env-var")
  expect_true(isTRUE(captured_warning$makevars_user_exists))
  expect_true(isTRUE(captured_warning$makevars_user_file_empty))
  expect_false(isTRUE(captured_warning$makevars_user_explicit_missing))
  expect_false(isTRUE(captured_warning$makevars_user_has_toolchain_override))
  expect_false(isTRUE(captured_warning$makevars_override_clears_stale_flibs))
  expect_identical(
    captured_warning$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    captured_warning$runtime_hint,
    "configured R_MAKEVARS_USER file is present but empty",
    fixed = TRUE
  )
  expect_false(grepl(
    "discovered user Makevars file",
    captured_warning$runtime_hint,
    fixed = TRUE
  ))
})

test_that("TC-9.4.26: non-empty explicit R_MAKEVARS_USER gets a dedicated no-override hint", {
  home_dir <- tempfile("wcb-home-nonempty-explicit-override-")
  makevars_user_path <- file.path(home_dir, "override.mk")
  dir.create(home_dir, recursive = TRUE)
  writeLines(
    c(
      "CFLAGS = -O2",
      "PKG_CPPFLAGS = -DUSE_STORY_WORKER_PROBE"
    ),
    makevars_user_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_identical(
    runtime_diagnostics$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(runtime_diagnostics$makevars_user_source, "env-var")
  expect_true(isTRUE(runtime_diagnostics$makevars_user_exists))
  expect_false(isTRUE(runtime_diagnostics$makevars_user_file_empty))
  expect_false(isTRUE(runtime_diagnostics$makevars_user_explicit_missing))
  expect_false(isTRUE(runtime_diagnostics$makevars_user_has_toolchain_override))
  expect_false(isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "configured R_MAKEVARS_USER file is non-empty",
    fixed = TRUE
  )
  expect_false(grepl(
    "discovered user Makevars file",
    runtime_diagnostics$runtime_hint,
    fixed = TRUE
  ))
})

test_that("TC-9.4.26: fwildclusterboot warnings surface non-empty explicit R_MAKEVARS_USER diagnostics", {
  home_dir <- tempfile("wcb-home-nonempty-explicit-override-warning-")
  makevars_user_path <- file.path(home_dir, "override.mk")
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  dir.create(home_dir, recursive = TRUE)
  writeLines(
    c(
      "CFLAGS = -O2",
      "PKG_CPPFLAGS = -DUSE_STORY_WORKER_WARNING_PROBE"
    ),
    makevars_user_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(
    captured_warning$makevars_user_path,
    normalizePath(makevars_user_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(captured_warning$makevars_user_source, "env-var")
  expect_true(isTRUE(captured_warning$makevars_user_exists))
  expect_false(isTRUE(captured_warning$makevars_user_file_empty))
  expect_false(isTRUE(captured_warning$makevars_user_explicit_missing))
  expect_false(isTRUE(captured_warning$makevars_user_has_toolchain_override))
  expect_false(isTRUE(captured_warning$makevars_override_clears_stale_flibs))
  expect_identical(
    captured_warning$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    captured_warning$runtime_hint,
    "configured R_MAKEVARS_USER file is non-empty",
    fixed = TRUE
  )
  expect_false(grepl(
    "discovered user Makevars file",
    captured_warning$runtime_hint,
    fixed = TRUE
  ))
})

test_that("TC-9.4.26: missing explicit R_MAKEVARS_USER is flagged as a missing override", {
  home_dir <- tempfile("wcb-home-missing-override-")
  makevars_dir <- file.path(home_dir, ".R")
  default_makevars_path <- file.path(makevars_dir, "Makevars")
  missing_makevars_path <- file.path(home_dir, "missing", "override.mk")
  dir.create(makevars_dir, recursive = TRUE)
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    default_makevars_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.setenv(R_MAKEVARS_USER = missing_makevars_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_identical(
    runtime_diagnostics$makevars_user_path,
    normalizePath(missing_makevars_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(runtime_diagnostics$makevars_user_source, "env-var")
  expect_false(isTRUE(runtime_diagnostics$makevars_user_exists))
  expect_true(isTRUE(runtime_diagnostics$makevars_user_explicit_missing))
  expect_false(isTRUE(runtime_diagnostics$makevars_override_clears_stale_flibs))
  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "configured R_MAKEVARS_USER path does not exist",
    fixed = TRUE
  )
  expect_match(
    runtime_diagnostics$runtime_hint,
    "unset R_MAKEVARS_USER",
    fixed = TRUE
  )
})

test_that("TC-9.4.26: fwildclusterboot warnings surface missing explicit override diagnostics", {
  home_dir <- tempfile("wcb-home-missing-override-warning-")
  makevars_dir <- file.path(home_dir, ".R")
  default_makevars_path <- file.path(makevars_dir, "Makevars")
  missing_makevars_path <- file.path(home_dir, "missing", "override.mk")
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  dir.create(makevars_dir, recursive = TRUE)
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    default_makevars_path,
    useBytes = TRUE
  )

  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_home)) {
      Sys.unsetenv("HOME")
    } else {
      Sys.setenv(HOME = old_home)
    }
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(home_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(HOME = home_dir)
  Sys.setenv(R_MAKEVARS_USER = missing_makevars_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(
    captured_warning$makevars_user_path,
    normalizePath(missing_makevars_path, winslash = "/", mustWork = FALSE)
  )
  expect_identical(captured_warning$makevars_user_source, "env-var")
  expect_false(isTRUE(captured_warning$makevars_user_exists))
  expect_true(isTRUE(captured_warning$makevars_user_explicit_missing))
  expect_false(isTRUE(captured_warning$makevars_override_clears_stale_flibs))
  expect_identical(
    captured_warning$blocker_boundary,
    "local-runtime-missing-fwildclusterboot-toolchain"
  )
  expect_match(
    captured_warning$runtime_hint,
    "configured R_MAKEVARS_USER path does not exist",
    fixed = TRUE
  )
  expect_match(
    captured_warning$runtime_hint,
    "unset R_MAKEVARS_USER",
    fixed = TRUE
  )
})

test_that("TC-9.4.26: fwildclusterboot warnings surface Makevars override diagnostics", {
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  runtime_diagnostics <- list(
    blocker_boundary = "d13-skip-backed-optional-backend-watchpoint",
    toolchain_mismatch_detected = TRUE,
    flibs_missing_paths = "/opt/gfortran/lib",
    homebrew_gfortran_candidates =
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib",
    makeconf_path = "/Library/Frameworks/R.framework/Resources/etc/Makeconf",
    makeconf_exists = TRUE,
    makeconf_flibs =
      "-L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
    makeconf_fc = "/opt/gfortran/bin/gfortran",
    makeconf_f77 = "/opt/gfortran/bin/gfortran",
    makeconf_uses_stale_opt_gfortran = TRUE,
    makevars_user_path = "/tmp/lwdid-makevars.mk",
    makevars_user_exists = TRUE,
    makevars_user_source = "env-var",
    makevars_user_file_empty = FALSE,
    makevars_user_has_toolchain_override = TRUE,
    makevars_user_flibs =
      "-L/opt/homebrew/lib/gcc/current -lgfortran -lemutls_w -lquadmath",
    makevars_user_fc = "/opt/homebrew/bin/gfortran",
    makevars_user_f77 = "/opt/homebrew/bin/gfortran",
    makevars_override_clears_stale_flibs = TRUE,
    runtime_hint = paste(
      "The configured user Makevars clears the stale local FLIBS floor,",
      "so the remaining fwildclusterboot unblocker should be tracked as",
      "package-level dependency/source-closure work."
    )
  )

  testthat::local_mocked_bindings(
    .wcb_collect_fwildclusterboot_runtime_diagnostics = function() runtime_diagnostics,
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(
    captured_warning$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )
  expect_identical(
    captured_warning$makevars_user_path,
    runtime_diagnostics$makevars_user_path
  )
  expect_identical(
    captured_warning$makevars_user_exists,
    runtime_diagnostics$makevars_user_exists
  )
  expect_identical(
    captured_warning$makevars_user_source,
    runtime_diagnostics$makevars_user_source
  )
  expect_identical(
    captured_warning$makevars_user_file_empty,
    runtime_diagnostics$makevars_user_file_empty
  )
  expect_identical(
    captured_warning$makevars_user_has_toolchain_override,
    runtime_diagnostics$makevars_user_has_toolchain_override
  )
  expect_identical(
    captured_warning$makevars_user_flibs,
    runtime_diagnostics$makevars_user_flibs
  )
  expect_identical(
    captured_warning$makevars_user_fc,
    runtime_diagnostics$makevars_user_fc
  )
  expect_identical(
    captured_warning$makevars_user_f77,
    runtime_diagnostics$makevars_user_f77
  )
  expect_identical(
    captured_warning$makevars_override_clears_stale_flibs,
    runtime_diagnostics$makevars_override_clears_stale_flibs
  )
})

test_that("TC-9.4.26: package-source-closure warnings surface dependency status payload", {
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  runtime_diagnostics <- list(
    blocker_boundary = "d13-skip-backed-optional-backend-watchpoint",
    toolchain_mismatch_detected = TRUE,
    flibs_missing_paths = "/opt/gfortran/lib",
    homebrew_gfortran_candidates =
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib",
    makeconf_path = "/Library/Frameworks/R.framework/Resources/etc/Makeconf",
    makeconf_exists = TRUE,
    makeconf_flibs =
      "-L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
    makeconf_fc = "/opt/gfortran/bin/gfortran",
    makeconf_f77 = "/opt/gfortran/bin/gfortran",
    makeconf_uses_stale_opt_gfortran = TRUE,
    makevars_user_path = "/tmp/lwdid-makevars.mk",
    makevars_user_exists = TRUE,
    makevars_user_source = "env-var",
    makevars_user_file_empty = FALSE,
    makevars_user_has_toolchain_override = TRUE,
    makevars_user_flibs =
      "-L/opt/homebrew/lib/gcc/current -lgfortran -lemutls_w -lquadmath",
    makevars_user_fc = "/opt/homebrew/bin/gfortran",
    makevars_user_f77 = "/opt/homebrew/bin/gfortran",
    makevars_override_clears_stale_flibs = TRUE,
    fwildclusterboot_available = FALSE,
    dependency_status = list(
      fwildclusterboot = FALSE,
      summclust = TRUE,
      JuliaConnectoR = TRUE,
      sitmo = FALSE,
      dqrng = FALSE
    ),
    missing_dependency_names = c("fwildclusterboot", "sitmo", "dqrng"),
    missing_dependency_count = 3L,
    direct_install_failure_node =
      "r-universe-index-visible-source-archive-403",
    direct_install_probe_provider = "r-universe",
    direct_install_provider_index_visible = TRUE,
    direct_install_source_archive_forbidden_detected = TRUE,
    direct_install_source_archive_status_code = 403L,
    github_direct_install_tested = FALSE,
    runtime_hint = paste(
      "The configured user Makevars clears the stale local FLIBS floor,",
      "so the remaining fwildclusterboot unblocker should be tracked as the",
      "r-universe-index-visible / r-universe-source-archive-403",
      "direct-install branch."
    )
  )

  testthat::local_mocked_bindings(
    .wcb_collect_fwildclusterboot_runtime_diagnostics = function() runtime_diagnostics,
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(
    captured_warning$fwildclusterboot_available,
    runtime_diagnostics$fwildclusterboot_available
  )
  expect_identical(
    captured_warning$dependency_status,
    runtime_diagnostics$dependency_status
  )
  expect_identical(
    captured_warning$missing_dependency_names,
    runtime_diagnostics$missing_dependency_names
  )
  expect_identical(
    captured_warning$missing_dependency_count,
    runtime_diagnostics$missing_dependency_count
  )
  expect_identical(
    captured_warning$direct_install_failure_node,
    runtime_diagnostics$direct_install_failure_node
  )
  expect_identical(
    captured_warning$direct_install_probe_provider,
    runtime_diagnostics$direct_install_probe_provider
  )
  expect_identical(
    captured_warning$direct_install_provider_index_visible,
    runtime_diagnostics$direct_install_provider_index_visible
  )
  expect_identical(
    captured_warning$direct_install_source_archive_forbidden_detected,
    runtime_diagnostics$direct_install_source_archive_forbidden_detected
  )
  expect_identical(
    captured_warning$direct_install_source_archive_status_code,
    runtime_diagnostics$direct_install_source_archive_status_code
  )
  expect_identical(
    captured_warning$github_direct_install_tested,
    runtime_diagnostics$github_direct_install_tested
  )
})

test_that("TC-9.4.26: actual package-source-closure warnings keep provider-scope status booleans machine-readable", {
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  makevars_user_path <- tempfile("wcb-makevars-provider-", fileext = ".mk")
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    makevars_user_path,
    useBytes = TRUE
  )
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(makevars_user_path)
  }, add = TRUE)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .fwildclusterboot_available = function() FALSE,
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_true(isTRUE(runtime_diagnostics$direct_install_provider_index_visible))
  expect_true(isTRUE(
    runtime_diagnostics$direct_install_source_archive_forbidden_detected
  ))
  expect_identical(
    runtime_diagnostics$direct_install_source_archive_status_code,
    403L
  )
  expect_identical(
    captured_warning$direct_install_provider_index_visible,
    runtime_diagnostics$direct_install_provider_index_visible
  )
  expect_identical(
    captured_warning$direct_install_source_archive_forbidden_detected,
    runtime_diagnostics$direct_install_source_archive_forbidden_detected
  )
  expect_identical(
    captured_warning$direct_install_source_archive_status_code,
    runtime_diagnostics$direct_install_source_archive_status_code
  )
})

test_that("TC-9.4.13: fwildclusterboot runtime-failure warning carries source error metadata", {
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  testthat::local_mocked_bindings(
    .fwildclusterboot_available = function() TRUE,
    .wcb_via_fwildclusterboot = function(...) {
      stop("simulated fwildclusterboot adapter failure", call. = FALSE)
    },
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(captured_warning$detail, "fwildclusterboot_runtime_failure")
  expect_identical(
    captured_warning$action_taken,
    "falling back to native bootstrap path"
  )
  expect_identical(
    captured_warning$source_error,
    "simulated fwildclusterboot adapter failure"
  )
  expect_identical(captured_warning$requested_n_bootstrap, 399L)
  expect_identical(captured_warning$actual_n_bootstrap, 399L)
  expect_false(isTRUE(captured_warning$full_enumeration))
  expect_identical(captured_warning$weight_type, "webb")
  expect_true(isTRUE(captured_warning$adapter_available))
  expect_identical(
    captured_warning$blocker_boundary,
    runtime_diagnostics$blocker_boundary
  )
  expect_identical(
    captured_warning$toolchain_mismatch_detected,
    runtime_diagnostics$toolchain_mismatch_detected
  )
  expect_identical(
    captured_warning$flibs_missing_paths,
    runtime_diagnostics$flibs_missing_paths
  )
  expect_identical(
    captured_warning$homebrew_gfortran_candidates,
    runtime_diagnostics$homebrew_gfortran_candidates
  )
  expect_identical(
    captured_warning$makeconf_path,
    runtime_diagnostics$makeconf_path
  )
  expect_identical(
    captured_warning$makeconf_exists,
    runtime_diagnostics$makeconf_exists
  )
  expect_identical(
    captured_warning$makeconf_flibs,
    runtime_diagnostics$makeconf_flibs
  )
  expect_identical(
    captured_warning$makeconf_fc,
    runtime_diagnostics$makeconf_fc
  )
  expect_identical(
    captured_warning$makeconf_f77,
    runtime_diagnostics$makeconf_f77
  )
  expect_identical(
    captured_warning$makeconf_uses_stale_opt_gfortran,
    runtime_diagnostics$makeconf_uses_stale_opt_gfortran
  )
  expect_identical(
    captured_warning$runtime_hint,
    runtime_diagnostics$runtime_hint
  )
})

test_that("TC-9.4.13: fwildclusterboot runtime-failure warnings keep the direct-install failure node machine-readable", {
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  runtime_diagnostics <- list(
    blocker_boundary = "d13-skip-backed-optional-backend-watchpoint",
    toolchain_mismatch_detected = TRUE,
    flibs_missing_paths = "/opt/gfortran/lib",
    homebrew_gfortran_candidates =
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib",
    makeconf_path = "/Library/Frameworks/R.framework/Resources/etc/Makeconf",
    makeconf_exists = TRUE,
    makeconf_flibs =
      "-L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
    makeconf_fc = "/opt/gfortran/bin/gfortran",
    makeconf_f77 = "/opt/gfortran/bin/gfortran",
    makeconf_uses_stale_opt_gfortran = TRUE,
    makevars_user_path = "/tmp/lwdid-makevars.mk",
    makevars_user_exists = TRUE,
    makevars_user_source = "env-var",
    makevars_user_file_empty = FALSE,
    makevars_user_has_toolchain_override = TRUE,
    makevars_user_flibs =
      "-L/opt/homebrew/lib/gcc/current -lgfortran -lemutls_w -lquadmath",
    makevars_user_fc = "/opt/homebrew/bin/gfortran",
    makevars_user_f77 = "/opt/homebrew/bin/gfortran",
    makevars_override_clears_stale_flibs = TRUE,
    fwildclusterboot_available = FALSE,
    dependency_status = list(
      fwildclusterboot = FALSE,
      summclust = TRUE,
      JuliaConnectoR = TRUE,
      sitmo = FALSE,
      dqrng = FALSE
    ),
    missing_dependency_names = c("fwildclusterboot", "sitmo", "dqrng"),
    missing_dependency_count = 3L,
    direct_install_failure_node =
      "r-universe-index-visible-source-archive-403",
    direct_install_probe_provider = "r-universe",
    direct_install_provider_index_visible = TRUE,
    direct_install_source_archive_forbidden_detected = TRUE,
    direct_install_source_archive_status_code = 403L,
    github_direct_install_tested = FALSE,
    runtime_hint = paste(
      "The configured user Makevars clears the stale local FLIBS floor,",
      "so the remaining fwildclusterboot unblocker should be tracked as the",
      "r-universe-index-visible / r-universe-source-archive-403",
      "direct-install branch."
    )
  )

  testthat::local_mocked_bindings(
    .wcb_collect_fwildclusterboot_runtime_diagnostics = function() runtime_diagnostics,
    .fwildclusterboot_available = function() TRUE,
    .wcb_via_fwildclusterboot = function(...) {
      stop("simulated fwildclusterboot adapter failure", call. = FALSE)
    },
    .package = "lwdid"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_identical(
    captured_warning$direct_install_failure_node,
    runtime_diagnostics$direct_install_failure_node
  )
  expect_identical(
    captured_warning$direct_install_probe_provider,
    runtime_diagnostics$direct_install_probe_provider
  )
  expect_identical(
    captured_warning$direct_install_provider_index_visible,
    runtime_diagnostics$direct_install_provider_index_visible
  )
  expect_identical(
    captured_warning$direct_install_source_archive_forbidden_detected,
    runtime_diagnostics$direct_install_source_archive_forbidden_detected
  )
  expect_identical(
    captured_warning$direct_install_source_archive_status_code,
    runtime_diagnostics$direct_install_source_archive_status_code
  )
  expect_true(isTRUE(captured_warning$github_direct_install_tested))
  expect_true(isTRUE(captured_warning$source_closure_cleared_for_github))
  expect_identical(
    captured_warning$blocker_boundary,
    "fwildclusterboot-github-runtime-failure"
  )
  expect_match(captured_warning$runtime_hint, "adapter runtime failure", fixed = TRUE)
})

test_that("TC-9.4.13: actual runtime-failure warnings keep provider-scope status booleans machine-readable", {
  fixture <- make_wcb_controls_fixture()
  captured_warning <- NULL
  makevars_user_path <- tempfile("wcb-runtime-provider-", fileext = ".mk")
  writeLines(
    c(
      "FLIBS = -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc -lgfortran -lemutls_w -lquadmath",
      "FC = /opt/homebrew/bin/gfortran",
      "F77 = /opt/homebrew/bin/gfortran"
    ),
    makevars_user_path,
    useBytes = TRUE
  )
  old_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  on.exit({
    if (is.na(old_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = old_makevars_user)
    }
    unlink(makevars_user_path)
  }, add = TRUE)
  Sys.setenv(R_MAKEVARS_USER = makevars_user_path)

  testthat::local_mocked_bindings(
    .wcb_get_flibs = function() {
      "-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"
    },
    .wcb_find_homebrew_gfortran_candidates = function() {
      "/opt/homebrew/Cellar/gcc/15.2.0/lib/gcc/current/libgfortran.5.dylib"
    },
    .wcb_get_makeconf_path = function() {
      "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
    },
    .wcb_read_makeconf_lines = function(path) {
      expect_identical(
        path,
        "/Library/Frameworks/R.framework/Resources/etc/Makeconf"
      )
      c(
        "FLIBS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0 -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath",
        "FC = /opt/gfortran/bin/gfortran",
        "F77 = /opt/gfortran/bin/gfortran"
      )
    },
    .fwildclusterboot_available = function() TRUE,
    .wcb_via_fwildclusterboot = function(...) {
      stop("'data' must be a data.frame, environment, or list", call. = FALSE)
    },
    .package = "lwdid"
  )

  runtime_diagnostics <- get(
    ".wcb_collect_fwildclusterboot_runtime_diagnostics",
    envir = asNamespace("lwdid"),
    inherits = FALSE
  )()

  expect_identical(
    runtime_diagnostics$blocker_boundary,
    "d13-skip-backed-optional-backend-watchpoint"
  )

  result <- withCallingHandlers(
    local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        controls = "x1",
        n_bootstrap = 399L,
        weight_type = "webb",
        alpha = 0.10,
        seed = 99L,
        impose_null = FALSE,
        use_fwildclusterboot = TRUE
      )
    ),
    lwdid_data = function(w) {
      captured_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_s3_class(captured_warning, "lwdid_data")
  expect_true(isTRUE(runtime_diagnostics$direct_install_provider_index_visible))
  expect_true(isTRUE(
    runtime_diagnostics$direct_install_source_archive_forbidden_detected
  ))
  expect_identical(
    runtime_diagnostics$direct_install_source_archive_status_code,
    403L
  )
  expect_identical(
    captured_warning$direct_install_provider_index_visible,
    runtime_diagnostics$direct_install_provider_index_visible
  )
  expect_identical(
    captured_warning$direct_install_source_archive_forbidden_detected,
    runtime_diagnostics$direct_install_source_archive_forbidden_detected
  )
  expect_identical(
    captured_warning$direct_install_source_archive_status_code,
    runtime_diagnostics$direct_install_source_archive_status_code
  )
  expect_identical(
    captured_warning$source_error,
    "'data' must be a data.frame, environment, or list"
  )
  expect_identical(
    captured_warning$blocker_boundary,
    "fwildclusterboot-github-runtime-failure"
  )
  expect_true(isTRUE(captured_warning$github_direct_install_tested))
  expect_true(isTRUE(captured_warning$source_closure_cleared_for_github))
  expect_match(captured_warning$runtime_hint, "adapter runtime failure", fixed = TRUE)
})

test_that("TC-9.4.24: .compute_bootstrap_pvalue_at_null() returns one at the observed ATT", {
  fixture <- make_wcb_controls_fixture()
  original_result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "with_controls"
    )
  )

  pvalue_at_att <- lwdid:::.compute_bootstrap_pvalue_at_null(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    null_value = original_result$att,
    att_original = original_result$att,
    se_original = original_result$original_se,
    n_bootstrap = 999L,
    weight_type = "rademacher",
    seed = 123L
  )

  expect_equal(pvalue_at_att, 1, tolerance = 1e-12)
})

test_that("TC-9.4.25: .compute_bootstrap_pvalue_at_null() matches the intercept-only replay away from the ATT", {
  fixture <- make_wcb_controls_fixture()
  original_result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "with_controls"
    )
  )
  null_value <- original_result$att + 10 * original_result$original_se
  adjusted_fixture <- fixture
  adjusted_fixture$Y <- fixture$Y - null_value * fixture$D
  replay_result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = adjusted_fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "intercept_only",
      use_fwildclusterboot = FALSE
    )
  )
  manual_expected <- mean(
    abs(replay_result$t_stats_bootstrap[is.finite(replay_result$t_stats_bootstrap)]) >=
      abs((original_result$att - null_value) / original_result$original_se)
  )

  helper_pvalue <- lwdid:::.compute_bootstrap_pvalue_at_null(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    null_value = null_value,
    att_original = original_result$att,
    se_original = original_result$original_se,
    n_bootstrap = 999L,
    weight_type = "rademacher",
    seed = 123L
  )

  expect_equal(helper_pvalue, manual_expected, tolerance = 1e-12)
  expect_lt(helper_pvalue, 0.1)
})

test_that("TC-9.4.15: public test inversion returns a finite CI that contains the ATT", {
  fixture <- utils::read.csv(
    "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e9_04_layer2_test_inversion_shared_controls_fixture.csv",
    stringsAsFactors = FALSE
  )

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_equal(result$ci_method, "test_inversion")
  expect_equal(result$weight_type, "rademacher")
  expect_equal(result$requested_n_bootstrap, 999L)
  expect_equal(result$actual_n_bootstrap, 999L)
  expect_false(isTRUE(result$full_enumeration))
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se_bootstrap))
  expect_true(is.finite(result$pvalue))
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
  expect_lte(result$ci_lower, result$att)
  expect_gte(result$ci_upper, result$att)
  expect_equal(result$att, -0.05445544554455495, tolerance = 1e-12)
  expect_equal(result$ci_lower, -0.9375382205743829, tolerance = 0.12)
  expect_equal(result$ci_upper, 0.7369054495575496, tolerance = 0.05)
  expect_equal(result$se_bootstrap, 0.42715399748263577, tolerance = 0.05)
  expect_gt(result$pvalue, 0.75)
  expect_lt(abs(result$pvalue - 0.8598598598598599), 0.08)
})

test_that("TC-9.4.16: public test inversion CI width stays in the same order as percentile_t", {
  fixture <- utils::read.csv(
    "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e9_04_layer2_test_inversion_shared_controls_fixture.csv",
    stringsAsFactors = FALSE
  )

  test_inversion <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L
    )
  )
  percentile_t <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = "x1",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      restricted_model = "intercept_only",
      use_fwildclusterboot = FALSE
    )
  )

  test_inversion_width <- test_inversion$ci_upper - test_inversion$ci_lower
  percentile_t_width <- percentile_t$ci_upper - percentile_t$ci_lower
  width_ratio <- test_inversion_width / percentile_t_width
  relative_width_gap <- abs(test_inversion_width - percentile_t_width) / abs(percentile_t_width)

  expect_s3_class(test_inversion, "lwdid_wcb_result")
  expect_s3_class(percentile_t, "lwdid_wcb_result")
  expect_identical(test_inversion$ci_method, "test_inversion")
  expect_identical(percentile_t$ci_method, "percentile_t")
  expect_identical(test_inversion$actual_n_bootstrap, 999L)
  expect_identical(percentile_t$actual_n_bootstrap, 32L)
  expect_false(isTRUE(test_inversion$full_enumeration))
  expect_true(isTRUE(percentile_t$full_enumeration))
  expect_true(is.finite(test_inversion_width))
  expect_true(is.finite(percentile_t_width))
  expect_gt(test_inversion_width, 0)
  expect_gt(percentile_t_width, 0)
  expect_gte(width_ratio, 0.5)
  expect_lte(width_ratio, 1.5)
  expect_lte(relative_width_gap, 0.5)
})

test_that("TC-9.4.18: null-effect WCB keeps most p-values above alpha on the story-local Monte Carlo DGP", {
  scenario <- make_wcb_story_local_decision_scenario(
    scenario_id = "story_local_null_effect_non_rejection",
    true_tau = 0,
    n_simulations = 100L,
    metric_name = "non_rejection_rate",
    decision_direction = "gt"
  )

  result <- lwdid:::.run_clustering_monte_carlo_wild_bootstrap_decision_scenario(
    scenario = scenario,
    shared_dgp = make_wcb_story_local_shared_dgp()
  )

  expect_identical(result$metric_name, "non_rejection_rate")
  expect_identical(result$requested_n_bootstrap, 199L)
  expect_identical(result$actual_n_bootstrap, 1024L)
  expect_identical(result$weight_type, "rademacher")
  expect_true(isTRUE(result$impose_null))
  expect_true(isTRUE(result$full_enumeration))
  expect_identical(result$ci_method, "percentile_t")
  expect_gt(result$metric_value, 0.88)
})

test_that("TC-9.4.19: story-local measured-blocker diagnostics stay executable and below the current threshold", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  expected_interval <- stats::binom.test(x = 16L, n = 50L)$conf.int

  expect_identical(diagnostics$case_id, "TC-9.4.19")
  expect_identical(
    diagnostics$exact_status,
    "story-local-power-gap-measured"
  )
  expect_identical(
    diagnostics$blocker_boundary,
    "story-local-hardening-power-gap-tc-9-4-19"
  )
  expect_identical(
    diagnostics$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(
    diagnostics$d13_runtime_classification,
    "d13-skip-backed-optional-backend-watchpoint"
  )
  expect_false(isTRUE(diagnostics$dedicated_suite_ready))
  expect_match(
    diagnostics$remaining_gap,
    "keep the dedicated-suite coverage landed",
    fixed = TRUE
  )
  expect_identical(diagnostics$scenario$metric_name, "power_under_alternative")
  expect_identical(diagnostics$scenario$decision_direction, "lt")
  expect_equal(diagnostics$scenario$true_tau, 2, tolerance = 0)
  expect_identical(diagnostics$scenario$G, 10L)
  expect_identical(diagnostics$scenario$obs_per_cluster, 20L)
  expect_identical(diagnostics$scenario$n_simulations, 50L)
  expect_equal(diagnostics$threshold, 0.70, tolerance = 0)

  expect_identical(diagnostics$result$metric_name, "power_under_alternative")
  expect_identical(diagnostics$result$n_simulations, 50L)
  expect_identical(diagnostics$result$n_attempts, 50L)
  expect_identical(diagnostics$result$rejection_count, 16L)
  expect_identical(diagnostics$result$requested_n_bootstrap, 199L)
  expect_identical(diagnostics$result$actual_n_bootstrap, 1024L)
  expect_true(isTRUE(diagnostics$result$impose_null))
  expect_true(isTRUE(diagnostics$result$full_enumeration))
  expect_identical(diagnostics$result$ci_method, "percentile_t")
  expect_type(diagnostics$result$replication_trace, "list")
  expect_identical(
    diagnostics$result$replication_trace$hit_count,
    diagnostics$result$rejection_count
  )
  expect_identical(
    diagnostics$result$replication_trace$miss_count,
    34L
  )
  expect_identical(
    length(diagnostics$result$replication_trace$hit_indices),
    diagnostics$result$rejection_count
  )
  expect_identical(
    length(diagnostics$result$replication_trace$miss_indices),
    diagnostics$result$replication_trace$miss_count
  )
  expect_identical(
    sort(c(
      diagnostics$result$replication_trace$hit_indices,
      diagnostics$result$replication_trace$miss_indices
    )),
    seq_len(diagnostics$result$n_simulations)
  )
  expect_equal(
    diagnostics$result$power_confidence_interval$lower,
    unname(expected_interval[[1L]]),
    tolerance = 1e-12
  )
  expect_equal(
    diagnostics$result$power_confidence_interval$upper,
    unname(expected_interval[[2L]]),
    tolerance = 1e-12
  )
  expect_identical(diagnostics$result$power_confidence_interval$conf_level, 0.95)
  expect_equal(diagnostics$result$metric_value, 0.32, tolerance = 0)
  expect_identical(diagnostics$required_rejection_count, 36L)
  expect_equal(diagnostics$required_power_to_clear, 0.72, tolerance = 0)
  expect_identical(diagnostics$rejection_count_shortfall, 20L)
  expect_false(isTRUE(diagnostics$threshold_passed))
  expect_equal(diagnostics$threshold_gap, -0.38, tolerance = 0)
})

test_that("story-local TC-9.4.19 diagnostics keep threshold prose in sync with custom hardening cutoffs", {
  below_threshold <- lwdid:::.story_local_wcb_power_diagnostics(threshold = 0.80)
  expect_identical(
    below_threshold$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_false(isTRUE(below_threshold$threshold_passed))
  expect_equal(below_threshold$threshold_gap, -0.48, tolerance = 0)
  expect_match(
    below_threshold$remaining_gap,
    "0.80 story-local threshold",
    fixed = TRUE
  )
  expect_match(
    below_threshold$remaining_gap,
    "dedicated-suite coverage landed",
    fixed = TRUE
  )

  above_threshold <- lwdid:::.story_local_wcb_power_diagnostics(threshold = 0.30)
  expect_true(isTRUE(above_threshold$threshold_passed))
  expect_true(isTRUE(above_threshold$dedicated_suite_ready))
  expect_length(above_threshold$global_blockers, 0L)
  expect_identical(
    above_threshold$story_live_blockers,
    "d13-skip-backed-optional-backend-watchpoint"
  )
  expect_identical(
    above_threshold$d13_runtime_classification,
    "d13-skip-backed-optional-backend-watchpoint"
  )
  expect_identical(above_threshold$required_rejection_count, 16L)
  expect_equal(above_threshold$required_power_to_clear, 0.32, tolerance = 0)
  expect_identical(above_threshold$rejection_count_shortfall, 0L)
  expect_equal(above_threshold$threshold_gap, 0.02, tolerance = 0)
  expect_match(
    above_threshold$remaining_gap,
    "0.30 story-local threshold",
    fixed = TRUE
  )
  expect_match(
    above_threshold$remaining_gap,
    "dedicated-suite coverage is already landed",
    fixed = TRUE
  )
})

test_that("story-local TC-9.4.19 diagnostics can rule out Mammen weights as an immediate power rescue", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(weight_type = "mammen")

  expect_identical(diagnostics$scenario$weight_type, "mammen")
  expect_identical(diagnostics$result$weight_type, "mammen")
  expect_identical(diagnostics$result$requested_n_bootstrap, 199L)
  expect_identical(diagnostics$result$actual_n_bootstrap, 199L)
  expect_false(isTRUE(diagnostics$result$full_enumeration))
  expect_identical(diagnostics$result$rejection_count, 3L)
  expect_equal(diagnostics$result$metric_value, 0.06, tolerance = 1e-12)
  expect_equal(diagnostics$threshold_gap, -0.64, tolerance = 1e-12)
  expect_identical(
    diagnostics$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
})

test_that("story-local TC-9.4.19 diagnostics can rule out Webb weights as the remaining simple power rescue", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(weight_type = "webb")

  expect_identical(diagnostics$scenario$weight_type, "webb")
  expect_identical(diagnostics$result$weight_type, "webb")
  expect_identical(diagnostics$result$requested_n_bootstrap, 199L)
  expect_identical(diagnostics$result$actual_n_bootstrap, 199L)
  expect_false(isTRUE(diagnostics$result$full_enumeration))
  expect_identical(diagnostics$result$rejection_count, 9L)
  expect_equal(diagnostics$result$metric_value, 0.18, tolerance = 1e-12)
  expect_equal(diagnostics$threshold_gap, -0.52, tolerance = 1e-12)
  expect_identical(
    diagnostics$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
})

test_that("story-local TC-9.4.19 weight-family screening compresses the simple distributions into one summary", {
  screening <- lwdid:::.story_local_wcb_weight_family_screening()

  expect_identical(
    screening$screened_weight_types,
    c("rademacher", "mammen", "webb")
  )
  expect_identical(
    vapply(screening$weight_family_results, `[[`, character(1L), "weight_type"),
    c("rademacher", "mammen", "webb")
  )
  expect_identical(
    vapply(screening$weight_family_results, `[[`, integer(1L), "rejection_count"),
    c(16L, 3L, 9L)
  )
  expect_equal(
    vapply(screening$weight_family_results, `[[`, numeric(1L), "metric_value"),
    c(0.32, 0.06, 0.18),
    tolerance = 1e-12
  )
  expect_identical(
    vapply(
      screening$weight_family_results,
      `[[`,
      integer(1L),
      "actual_n_bootstrap"
    ),
    c(1024L, 199L, 199L)
  )
  expect_identical(
    screening$comparison$ranking_best_to_worst,
    c("rademacher", "webb", "mammen")
  )
  expect_identical(screening$comparison$best_weight_type, "rademacher")
  expect_identical(
    screening$comparison$best_non_rademacher_weight_type,
    "webb"
  )
})

test_that("story-local TC-9.4.19 weight-family screening makes the shortcut verdict machine-readable", {
  screening <- lwdid:::.story_local_wcb_weight_family_screening()

  expect_identical(
    screening$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_true(isTRUE(screening$summary$canonical_rademacher_remains_best))
  expect_true(isTRUE(screening$summary$all_weight_families_fail_threshold))
  expect_true(isTRUE(screening$summary$non_rademacher_shortcut_eliminated))
  expect_true(isTRUE(screening$summary$full_enumeration_only_for_rademacher))
  expect_equal(
    screening$comparison$best_non_rademacher_minus_rademacher,
    -0.14,
    tolerance = 1e-12
  )
  expect_equal(
    screening$comparison$webb_minus_mammen,
    0.12,
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 diagnostics can replay the impose_null FALSE decomposition regime", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(impose_null = FALSE)
  expected_interval <- stats::binom.test(x = 0L, n = 50L)$conf.int

  expect_type(diagnostics$result$replication_trace, "list")
  expect_false(isTRUE(diagnostics$result$impose_null))
  expect_identical(diagnostics$result$ci_method, "percentile")
  expect_identical(diagnostics$result$rejection_count, 0L)
  expect_identical(diagnostics$result$replication_trace$hit_count, 0L)
  expect_identical(diagnostics$result$replication_trace$miss_count, 50L)
  expect_length(diagnostics$result$replication_trace$hit_indices, 0L)
  expect_identical(
    diagnostics$result$replication_trace$miss_indices,
    seq_len(diagnostics$result$n_simulations)
  )
  expect_equal(diagnostics$result$metric_value, 0, tolerance = 0)
  expect_equal(
    diagnostics$result$power_confidence_interval$lower,
    unname(expected_interval[[1L]]),
    tolerance = 1e-12
  )
  expect_equal(
    diagnostics$result$power_confidence_interval$upper,
    unname(expected_interval[[2L]]),
    tolerance = 1e-12
  )
  expect_identical(diagnostics$result$power_confidence_interval$conf_level, 0.95)
  expect_false(isTRUE(diagnostics$threshold_passed))
  expect_equal(diagnostics$threshold_gap, -0.70, tolerance = 0)
})

test_that("story-local TC-9.4.19 diagnostics expose replication-level treated-cluster mix", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  replications <- diagnostics$result$replication_trace$replications
  treated_cluster_profile <-
    diagnostics$result$replication_trace$treated_cluster_profile

  expect_type(replications, "list")
  expect_length(replications, diagnostics$result$n_simulations)
  expect_identical(
    vapply(replications, `[[`, integer(1L), "replication_id"),
    seq_len(diagnostics$result$n_simulations)
  )
  expect_true(all(!is.na(vapply(replications, `[[`, logical(1L), "reject"))))
  expect_true(all(vapply(replications, function(row) is.numeric(row$pvalue), logical(1L))))
  expect_true(all(vapply(replications, function(row) is.numeric(row$att), logical(1L))))
  expect_true(all(vapply(replications, function(row) is.numeric(row$original_se), logical(1L))))
  expect_true(all(vapply(replications, function(row) is.numeric(row$t_stat_original), logical(1L))))
  expect_true(all(vapply(replications, function(row) is.numeric(row$treated_cluster_count), logical(1L))))
  expect_true(
    all(
      vapply(replications, function(row) row$treated_cluster_count >= 1L, logical(1L))
    )
  )
  expect_true(
    all(
      vapply(replications, function(row) row$treated_cluster_count <= diagnostics$scenario$G - 1L, logical(1L))
    )
  )
  expect_identical(
    sum(vapply(replications, `[[`, logical(1L), "reject")),
    diagnostics$result$rejection_count
  )

  expect_type(treated_cluster_profile, "list")
  expect_gt(length(treated_cluster_profile), 0L)
  expect_identical(
    vapply(treated_cluster_profile, `[[`, integer(1L), "treated_cluster_count"),
    sort(vapply(treated_cluster_profile, `[[`, integer(1L), "treated_cluster_count"))
  )
  expect_identical(
    sum(vapply(treated_cluster_profile, `[[`, integer(1L), "n_replications")),
    diagnostics$result$n_simulations
  )
  expect_identical(
    sum(vapply(treated_cluster_profile, `[[`, integer(1L), "rejection_count")),
    diagnostics$result$rejection_count
  )
  expect_equal(
    sum(vapply(
      treated_cluster_profile,
      function(bucket) bucket$rejection_rate * bucket$n_replications,
      numeric(1L)
    )),
    diagnostics$result$rejection_count,
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 diagnostics expose paired standard-t rescue cases", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  replications <- diagnostics$result$replication_trace$replications
  counterfactual_summary <-
    diagnostics$result$replication_trace$standard_t_counterfactual_summary

  expect_true(all(vapply(
    replications,
    function(row) is.logical(row$standard_t_reject),
    logical(1L)
  )))
  expect_true(all(vapply(
    replications,
    function(row) is.numeric(row$standard_t_pvalue),
    logical(1L)
  )))
  expect_identical(
    sum(vapply(replications, `[[`, logical(1L), "standard_t_reject")),
    21L
  )

  expect_type(counterfactual_summary, "list")
  expect_identical(counterfactual_summary$standard_t_rejection_count, 21L)
  expect_identical(counterfactual_summary$wild_bootstrap_rejection_count, 16L)
  expect_identical(counterfactual_summary$shared_rejection_count, 16L)
  expect_identical(counterfactual_summary$standard_only_rejection_count, 5L)
  expect_identical(counterfactual_summary$wcb_only_rejection_count, 0L)
  expect_identical(counterfactual_summary$paired_rejection_gap, 5L)
  expect_equal(counterfactual_summary$counterfactual_power, 0.42, tolerance = 1e-12)
  expect_identical(
    counterfactual_summary$required_rejection_count_to_clear_threshold,
    36L
  )
  expect_equal(
    counterfactual_summary$required_power_to_clear_threshold,
    0.72,
    tolerance = 1e-12
  )
  expect_identical(counterfactual_summary$residual_rejection_shortfall, 15L)
  expect_false(isTRUE(counterfactual_summary$counterfactual_clears_threshold))
  expect_identical(
    as.numeric(unname(unlist(counterfactual_summary$near_threshold_window))),
    c(0.04, 0.06)
  )
  expect_identical(counterfactual_summary$standard_only_near_threshold_count, 2L)
  expect_identical(
    counterfactual_summary$standard_only_non_near_threshold_count,
    3L
  )
  expect_identical(
    counterfactual_summary$standard_only_near_threshold_replay_seeds,
    c(42L, 80L)
  )
  expect_identical(
    counterfactual_summary$standard_only_non_near_threshold_replay_seeds,
    c(60L, 74L, 91L)
  )
  expect_identical(
    counterfactual_summary$standard_only_replay_seeds,
    c(42L, 60L, 74L, 80L, 91L)
  )

  standard_only_bucket_ranking <- counterfactual_summary$standard_only_bucket_ranking
  expect_type(standard_only_bucket_ranking, "list")
  expect_identical(
    vapply(standard_only_bucket_ranking, `[[`, integer(1L), "treated_cluster_count"),
    c(2L, 7L, 8L, 9L)
  )
  expect_identical(
    vapply(standard_only_bucket_ranking, `[[`, integer(1L), "count"),
    c(2L, 1L, 1L, 1L)
  )

  standard_only_cases <- counterfactual_summary$standard_only_cases
  expect_type(standard_only_cases, "list")
  expect_identical(
    vapply(standard_only_cases, `[[`, integer(1L), "attempt_id"),
    c(1L, 19L, 33L, 39L, 50L)
  )
  expect_identical(
    vapply(standard_only_cases, `[[`, integer(1L), "replay_seed"),
    c(42L, 60L, 74L, 80L, 91L)
  )
  expect_identical(
    vapply(standard_only_cases, `[[`, integer(1L), "treated_cluster_count"),
    c(8L, 2L, 9L, 7L, 2L)
  )
  expect_equal(
    vapply(standard_only_cases, `[[`, numeric(1L), "wcb_pvalue"),
    c(0.0576171875, 0.0947265625, 0.2666015625, 0.05859375, 0.1181640625),
    tolerance = 1e-12
  )
  expect_lt(
    max(abs(
      vapply(standard_only_cases, `[[`, numeric(1L), "standard_t_pvalue") - c(
        0.00319035679022905,
        0.0103470056514641,
        0.000116172065856869,
        0.0487317415562678,
        0.00925579091359756
      )
    )),
    1e-8
  )
})

test_that("story-local TC-9.4.19 diagnostics compress threshold-split rescue summary fields", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  counterfactual_summary <-
    diagnostics$result$replication_trace$standard_t_counterfactual_summary

  expect_equal(
    counterfactual_summary$standard_only_near_threshold_share,
    0.4,
    tolerance = 1e-12
  )
  expect_equal(
    counterfactual_summary$standard_only_non_near_threshold_share,
    0.6,
    tolerance = 1e-12
  )
  expect_equal(
    counterfactual_summary$standard_only_near_threshold_replication_share,
    0.04,
    tolerance = 1e-12
  )
  expect_equal(
    counterfactual_summary$standard_only_non_near_threshold_replication_share,
    0.06,
    tolerance = 1e-12
  )
  expect_identical(
    counterfactual_summary$closest_standard_only_seed_to_threshold,
    42L
  )
  expect_equal(
    counterfactual_summary$closest_standard_only_gap_to_threshold,
    0.0076171875,
    tolerance = 1e-12
  )
  expect_identical(
    counterfactual_summary$closest_non_near_standard_only_seed,
    60L
  )
  expect_equal(
    counterfactual_summary$closest_non_near_standard_only_gap_to_window_upper,
    0.0347265625,
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 diagnostics expose threshold-split rescue bucket rankings", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  counterfactual_summary <-
    diagnostics$result$replication_trace$standard_t_counterfactual_summary

  near_threshold_bucket_ranking <-
    counterfactual_summary$standard_only_near_threshold_bucket_ranking
  expect_type(near_threshold_bucket_ranking, "list")
  expect_identical(
    vapply(
      near_threshold_bucket_ranking,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ),
    c(7L, 8L)
  )
  expect_identical(
    vapply(near_threshold_bucket_ranking, `[[`, integer(1L), "count"),
    c(1L, 1L)
  )

  non_near_threshold_bucket_ranking <-
    counterfactual_summary$standard_only_non_near_threshold_bucket_ranking
  expect_type(non_near_threshold_bucket_ranking, "list")
  expect_identical(
    vapply(
      non_near_threshold_bucket_ranking,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ),
    c(2L, 9L)
  )
  expect_identical(
    vapply(non_near_threshold_bucket_ranking, `[[`, integer(1L), "count"),
    c(2L, 1L)
  )
})

test_that("story-local TC-9.4.19 standard-only rescue frontier orders cases by smallest WCB gap", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  frontier <-
    diagnostics$result$replication_trace$standard_t_counterfactual_summary$standard_only_repair_frontier

  expect_type(frontier, "list")
  expect_length(frontier, 5L)
  expect_identical(
    vapply(frontier, `[[`, integer(1L), "case_rank"),
    c(1L, 2L, 3L, 4L, 5L)
  )
  expect_identical(
    vapply(frontier, `[[`, integer(1L), "replay_seed"),
    c(42L, 80L, 60L, 91L, 74L)
  )
  expect_identical(
    vapply(frontier, `[[`, integer(1L), "replay_seed_offset"),
    c(0L, 38L, 18L, 49L, 32L)
  )
  expect_identical(
    vapply(frontier, `[[`, integer(1L), "attempt_id"),
    c(1L, 39L, 19L, 50L, 33L)
  )
  expect_identical(
    vapply(frontier, `[[`, integer(1L), "treated_cluster_count"),
    c(8L, 7L, 2L, 2L, 9L)
  )
  expect_equal(
    vapply(frontier, `[[`, numeric(1L), "wcb_pvalue_gap"),
    c(
      0.0076171875,
      0.00859375,
      0.0447265625,
      0.0681640625,
      0.2166015625
    ),
    tolerance = 1e-12
  )
  frontier_wcb_pvalues <- vapply(frontier, `[[`, numeric(1L), "wcb_pvalue")
  frontier_standard_t_pvalues <- vapply(
    frontier,
    `[[`,
    numeric(1L),
    "standard_t_pvalue"
  )
  expect_equal(
    vapply(frontier, `[[`, numeric(1L), "paired_pvalue_gap"),
    frontier_wcb_pvalues - frontier_standard_t_pvalues,
    tolerance = 1e-10
  )
})

test_that("story-local TC-9.4.19 standard-only rescue frontier still cannot clear threshold", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  counterfactual_summary <-
    diagnostics$result$replication_trace$standard_t_counterfactual_summary
  frontier <- counterfactual_summary$standard_only_repair_frontier

  expect_identical(
    vapply(frontier, `[[`, integer(1L), "cumulative_rescue_count"),
    c(1L, 2L, 3L, 4L, 5L)
  )
  expect_identical(
    vapply(frontier, `[[`, integer(1L), "counterfactual_rejection_count"),
    c(17L, 18L, 19L, 20L, 21L)
  )
  expect_equal(
    vapply(frontier, `[[`, numeric(1L), "counterfactual_power"),
    c(0.34, 0.36, 0.38, 0.40, 0.42),
    tolerance = 1e-12
  )
  expect_identical(
    vapply(frontier, `[[`, integer(1L), "residual_shortfall"),
    c(19L, 18L, 17L, 16L, 15L)
  )
  expect_identical(
    vapply(frontier, `[[`, logical(1L), "clears_threshold"),
    c(FALSE, FALSE, FALSE, FALSE, FALSE)
  )
  expect_equal(
    vapply(frontier, `[[`, numeric(1L), "cumulative_share_of_standard_only_cases"),
    c(0.2, 0.4, 0.6, 0.8, 1.0),
    tolerance = 1e-12
  )
  expect_null(counterfactual_summary$minimum_case_rank_to_clear_threshold)
  expect_null(counterfactual_summary$minimum_case_replay_seed_combo_to_clear_threshold)
  expect_null(counterfactual_summary$minimum_case_combo_share_of_standard_only_cases)
})

test_that("story-local TC-9.4.19 diagnostics expose bucket-level strength summaries", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  treated_cluster_profile <-
    diagnostics$result$replication_trace$treated_cluster_profile
  profile_summary <-
    diagnostics$result$replication_trace$treated_cluster_profile_summary
  profile_by_count <- stats::setNames(
    treated_cluster_profile,
    vapply(
      treated_cluster_profile,
      function(bucket) as.character(bucket$treated_cluster_count),
      character(1L)
    )
  )

  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$avg_abs_att),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$median_abs_att),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$avg_original_se),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$median_original_se),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$avg_pvalue),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$median_pvalue),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$avg_abs_t_stat),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$max_abs_t_stat),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$avg_abs_att_to_se_ratio),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$median_abs_att_to_se_ratio),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$min_abs_att_to_se_ratio),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.numeric(bucket$max_abs_att_to_se_ratio),
    logical(1L)
  )))
  expect_true(all(vapply(
    treated_cluster_profile,
    function(bucket) is.integer(bucket$near_threshold_count),
    logical(1L)
  )))

  expect_identical(profile_by_count[["5"]]$n_replications, 16L)
  expect_identical(profile_by_count[["5"]]$rejection_count, 3L)
  expect_lt(profile_by_count[["5"]]$avg_abs_att, profile_by_count[["3"]]$avg_abs_att)
  expect_gt(
    profile_by_count[["5"]]$avg_original_se,
    profile_by_count[["3"]]$avg_original_se
  )
  expect_gt(profile_by_count[["5"]]$avg_pvalue, profile_by_count[["3"]]$avg_pvalue)
  expect_lt(
    profile_by_count[["5"]]$avg_abs_t_stat,
    profile_by_count[["3"]]$avg_abs_t_stat
  )
  expect_equal(
    profile_by_count[["3"]]$avg_abs_att_to_se_ratio,
    2.413654960182192,
    tolerance = 1e-12
  )
  expect_equal(
    profile_by_count[["3"]]$median_abs_att_to_se_ratio,
    2.641230268234367,
    tolerance = 1e-12
  )
  expect_equal(
    profile_by_count[["3"]]$min_abs_att_to_se_ratio,
    0.719632578921078,
    tolerance = 1e-12
  )
  expect_equal(
    profile_by_count[["3"]]$max_abs_att_to_se_ratio,
    3.986240073811386,
    tolerance = 1e-12
  )
  expect_equal(
    profile_by_count[["5"]]$avg_abs_att_to_se_ratio,
    1.806857237921066,
    tolerance = 1e-12
  )
  expect_equal(
    profile_by_count[["5"]]$median_abs_att_to_se_ratio,
    1.500194139279165,
    tolerance = 1e-12
  )
  expect_equal(
    profile_by_count[["5"]]$min_abs_att_to_se_ratio,
    0.420512965833072,
    tolerance = 1e-12
  )
  expect_equal(
    profile_by_count[["5"]]$max_abs_att_to_se_ratio,
    4.498817868131472,
    tolerance = 1e-12
  )
  expect_lt(
    profile_by_count[["5"]]$avg_abs_att_to_se_ratio,
    profile_by_count[["3"]]$avg_abs_att_to_se_ratio
  )
  expect_lt(
    profile_by_count[["5"]]$median_abs_att_to_se_ratio,
    profile_by_count[["3"]]$median_abs_att_to_se_ratio
  )
  expect_identical(profile_by_count[["5"]]$near_threshold_count, 1L)

  expect_type(profile_summary, "list")
  expect_identical(profile_summary$modal_bucket$treated_cluster_count, 5L)
  expect_identical(profile_summary$modal_bucket$n_replications, 16L)
  expect_identical(profile_summary$modal_bucket$rejection_count, 3L)
  expect_equal(profile_summary$modal_bucket$rejection_rate, 3 / 16, tolerance = 1e-12)

  expect_identical(
    unlist(profile_summary$high_rejection_buckets$treated_cluster_counts),
    c(3L, 6L)
  )
  expect_identical(profile_summary$high_rejection_buckets$total_replications, 16L)
  expect_identical(profile_summary$high_rejection_buckets$total_rejections, 9L)
  expect_equal(
    profile_summary$high_rejection_buckets$share_of_rejections,
    9 / 16,
    tolerance = 1e-12
  )
  expect_equal(
    profile_summary$high_rejection_buckets$share_of_replications,
    16 / 50,
    tolerance = 1e-12
  )

  expect_identical(
    unlist(profile_summary$zero_rejection_buckets$treated_cluster_counts),
    c(2L, 8L, 9L)
  )
  expect_identical(profile_summary$zero_rejection_buckets$total_replications, 6L)

  expect_identical(
    unlist(profile_summary$near_threshold_buckets$treated_cluster_counts),
    c(4L, 5L, 6L, 7L, 8L)
  )
  expect_identical(
    profile_summary$near_threshold_buckets$total_near_threshold_cases,
    5L
  )

  contrast <- profile_summary$modal_bucket_signal_ratio_contrast
  expect_type(contrast, "list")
  expect_identical(contrast$modal_bucket_treated_cluster_count, 5L)
  expect_identical(contrast$comparison_bucket_treated_cluster_count, 3L)
  expect_true(isTRUE(contrast$comparison_bucket_is_high_rejection))
  expect_equal(
    contrast$modal_avg_abs_att_to_se_ratio,
    1.806857237921066,
    tolerance = 1e-12
  )
  expect_equal(
    contrast$comparison_avg_abs_att_to_se_ratio,
    2.413654960182192,
    tolerance = 1e-12
  )
  expect_equal(
    contrast$modal_median_abs_att_to_se_ratio,
    1.500194139279165,
    tolerance = 1e-12
  )
  expect_equal(
    contrast$comparison_median_abs_att_to_se_ratio,
    2.641230268234367,
    tolerance = 1e-12
  )
  expect_lt(contrast$avg_ratio_gap_modal_minus_comparison, 0)
  expect_lt(contrast$median_ratio_gap_modal_minus_comparison, 0)
  expect_true(isTRUE(contrast$modal_bucket_weaker_on_avg_ratio))
  expect_true(isTRUE(contrast$modal_bucket_weaker_on_median_ratio))
})

test_that("story-local TC-9.4.19 replication-mix helper freezes the modal-versus-anchor summary", {
  replication_mix <- lwdid:::.story_local_wcb_replication_mix_summary()

  expect_identical(replication_mix$case_id, "TC-9.4.19")
  expect_identical(
    replication_mix$exact_status,
    "story-local-power-gap-replication-mix-bridge-landed"
  )
  expect_identical(
    replication_mix$numeric_status,
    "treated-cluster-mix-diagnostics-frozen"
  )
  expect_identical(
    replication_mix$blocker_boundary,
    "story-local-hardening-power-gap-tc-9-4-19"
  )
  expect_identical(
    replication_mix$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_equal(replication_mix$threshold, 0.70, tolerance = 0)
  expect_identical(replication_mix$summary$total_replications, 50L)
  expect_identical(replication_mix$summary$total_rejections, 16L)
  expect_identical(replication_mix$summary$threshold_shortfall, 20L)
  expect_identical(
    replication_mix$summary$modal_bucket$treated_cluster_count,
    5L
  )
  expect_identical(
    replication_mix$summary$signal_ratio_anchor_bucket$treated_cluster_count,
    3L
  )
  expect_equal(
    replication_mix$summary$signal_ratio_advantage_over_modal$avg_abs_att_to_se_ratio,
    0.606797722261126,
    tolerance = 1e-12
  )
  expect_equal(
    replication_mix$summary$signal_ratio_advantage_over_modal$median_abs_att_to_se_ratio,
    1.141036128955202,
    tolerance = 1e-12
  )
  expect_identical(
    replication_mix$summary$high_rejection_buckets$treated_cluster_counts,
    c(3L, 6L)
  )
  expect_equal(
    replication_mix$summary$high_rejection_buckets$share_of_rejections,
    9 / 16,
    tolerance = 1e-12
  )
  expect_equal(
    replication_mix$summary$high_rejection_buckets$share_of_replications,
    16 / 50,
    tolerance = 1e-12
  )
  expect_identical(
    replication_mix$summary$zero_rejection_buckets$treated_cluster_counts,
    c(2L, 8L, 9L)
  )
  expect_identical(
    replication_mix$summary$near_threshold_buckets$treated_cluster_counts,
    c(4L, 5L, 6L, 7L, 8L)
  )
  expect_identical(
    replication_mix$summary$near_threshold_buckets$total_near_threshold_cases,
    5L
  )
})

test_that("story-local TC-9.4.19 signal-ratio bucket helper exposes modal and outlier rankings", {
  signal_ratio_profile <- lwdid:::.story_local_wcb_signal_ratio_bucket_profile()

  expect_identical(signal_ratio_profile$case_id, "TC-9.4.19")
  expect_identical(
    signal_ratio_profile$exact_status,
    "story-local-signal-ratio-bucket-profile-frozen"
  )
  expect_identical(
    signal_ratio_profile$numeric_status,
    "signal-ratio-bucket-profile-machine-readable"
  )
  expect_identical(
    signal_ratio_profile$blocker_boundary,
    "story-local-hardening-power-gap-tc-9-4-19"
  )
  expect_identical(
    signal_ratio_profile$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(
    signal_ratio_profile$consumer_summary_source,
    "treated_cluster_profile_summary"
  )
  expect_identical(
    signal_ratio_profile$profile_source,
    "treated_cluster_profile"
  )
  expect_identical(
    vapply(
      signal_ratio_profile$signal_ratio_ranking,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ),
    c(9L, 8L, 2L, 3L, 4L, 6L, 5L, 7L)
  )
  expect_equal(
    vapply(
      signal_ratio_profile$signal_ratio_ranking[1:3],
      `[[`,
      numeric(1L),
      "avg_abs_att_to_se_ratio"
    ),
    c(6.464390434806528, 2.928663161784076, 2.61096192545123),
    tolerance = 1e-12
  )
  expect_identical(
    signal_ratio_profile$summary$modal_bucket_rank_by_avg_ratio,
    7L
  )
  expect_identical(
    signal_ratio_profile$summary$comparison_bucket_rank_by_avg_ratio,
    4L
  )
  expect_identical(
    signal_ratio_profile$summary$strongest_bucket_overall$treated_cluster_count,
    9L
  )
  expect_identical(
    signal_ratio_profile$summary$strongest_high_rejection_bucket$treated_cluster_count,
    3L
  )
  expect_identical(
    signal_ratio_profile$summary$strongest_zero_rejection_bucket$treated_cluster_count,
    9L
  )
  expect_equal(
    signal_ratio_profile$summary$strongest_zero_rejection_bucket$avg_abs_att_to_se_ratio,
    6.464390434806528,
    tolerance = 1e-12
  )
  expect_equal(
    signal_ratio_profile$summary$strongest_zero_rejection_bucket$avg_pvalue,
    0.2666015625,
    tolerance = 0
  )
  expect_equal(
    signal_ratio_profile$summary$modal_vs_comparison_avg_ratio_gap,
    -0.606797722261126,
    tolerance = 1e-12
  )
  expect_equal(
    signal_ratio_profile$summary$modal_vs_comparison_median_ratio_gap,
    -1.141036128955202,
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 signal-ratio bucket helper crosswalks miss mass against ratio order", {
  signal_ratio_profile <- lwdid:::.story_local_wcb_signal_ratio_bucket_profile()

  expect_identical(
    vapply(
      signal_ratio_profile$signal_ratio_ranking,
      `[[`,
      integer(1L),
      "non_near_threshold_miss_count"
    ),
    c(1L, 1L, 3L, 2L, 4L, 5L, 12L, 3L)
  )
  expect_equal(
    vapply(
      signal_ratio_profile$signal_ratio_ranking,
      `[[`,
      numeric(1L),
      "share_of_non_near_threshold_misses"
    ),
    c(1, 1, 3, 2, 4, 5, 12, 3) / 31,
    tolerance = 1e-12
  )
  expect_identical(
    vapply(
      signal_ratio_profile$signal_ratio_ranking,
      `[[`,
      integer(1L),
      "non_near_threshold_miss_rank"
    ),
    c(8L, 7L, 4L, 6L, 3L, 2L, 1L, 5L)
  )
  expect_identical(
    signal_ratio_profile$summary$dominant_non_near_threshold_miss_bucket$treated_cluster_count,
    5L
  )
  expect_identical(
    signal_ratio_profile$summary$dominant_non_near_threshold_miss_bucket_rank_by_avg_ratio,
    7L
  )
  expect_equal(
    signal_ratio_profile$summary$dominant_non_near_threshold_miss_share,
    12 / 31,
    tolerance = 1e-12
  )
  expect_identical(
    signal_ratio_profile$summary$strongest_bucket_overall_non_near_threshold_miss_rank,
    8L
  )
  expect_equal(
    signal_ratio_profile$summary$strongest_bucket_overall_non_near_threshold_miss_share,
    1 / 31,
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 signal-ratio repair crosswalk isolates threshold-clearing combo", {
  crosswalk <- lwdid:::.story_local_wcb_signal_ratio_repair_crosswalk()

  expect_identical(crosswalk$case_id, "TC-9.4.19")
  expect_identical(
    crosswalk$exact_status,
    "story-local-signal-ratio-repair-crosswalk-frozen"
  )
  expect_identical(
    crosswalk$numeric_status,
    "signal-ratio-repair-crosswalk-machine-readable"
  )
  expect_identical(
    crosswalk$blocker_boundary,
    "story-local-hardening-power-gap-tc-9-4-19"
  )
  expect_identical(
    crosswalk$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(
    crosswalk$consumer_summary_source,
    "story_local_wcb_signal_ratio_bucket_profile()+bucket_repair_frontier"
  )

  strongest_ratio_bucket <- crosswalk$summary$strongest_ratio_bucket
  expect_identical(strongest_ratio_bucket$treated_cluster_count, 9L)
  expect_identical(strongest_ratio_bucket$signal_ratio_rank, 1L)
  expect_identical(strongest_ratio_bucket$repair_frontier_rank, 8L)
  expect_false(strongest_ratio_bucket$in_minimum_threshold_clearing_combo)
  expect_equal(
    strongest_ratio_bucket$share_of_non_near_threshold_misses,
    1 / 31,
    tolerance = 1e-12
  )

  minimum_combo <- crosswalk$summary$minimum_threshold_clearing_combo
  expect_type(minimum_combo, "list")
  expect_identical(
    vapply(minimum_combo, `[[`, integer(1L), "treated_cluster_count"),
    c(5L, 6L, 4L)
  )
  expect_identical(
    vapply(minimum_combo, `[[`, integer(1L), "signal_ratio_rank"),
    c(7L, 6L, 5L)
  )
  expect_identical(
    vapply(minimum_combo, `[[`, integer(1L), "repair_frontier_rank"),
    c(1L, 2L, 3L)
  )
  expect_identical(
    vapply(minimum_combo, `[[`, logical(1L), "in_minimum_threshold_clearing_combo"),
    c(TRUE, TRUE, TRUE)
  )
  expect_equal(
    vapply(minimum_combo, `[[`, numeric(1L), "counterfactual_power"),
    c(0.56, 0.66, 0.74),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(minimum_combo, `[[`, numeric(1L), "cumulative_share_of_non_near_threshold_misses"),
    c(12 / 31, 17 / 31, 21 / 31),
    tolerance = 1e-12
  )

  expect_identical(
    crosswalk$summary$first_threshold_clearing_bucket$treated_cluster_count,
    4L
  )
  expect_identical(
    crosswalk$summary$first_threshold_clearing_bucket$signal_ratio_rank,
    5L
  )
  expect_true(isTRUE(crosswalk$summary$strongest_ratio_bucket_excluded_from_repair_combo))
  expect_true(isTRUE(crosswalk$summary$threshold_clearing_combo_uses_mid_ratio_buckets))
})

test_that("story-local TC-9.4.19 power repair ladder orders the frozen counterfactual rungs", {
  ladder <- lwdid:::.story_local_wcb_power_repair_ladder()

  expect_identical(ladder$case_id, "TC-9.4.19")
  expect_identical(
    ladder$exact_status,
    "story-local-power-repair-ladder-frozen"
  )
  expect_identical(
    ladder$numeric_status,
    "power-repair-ladder-machine-readable"
  )
  expect_identical(
    ladder$blocker_boundary,
    "story-local-hardening-power-gap-tc-9-4-19"
  )
  expect_identical(
    ladder$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(
    ladder$consumer_summary_source,
    paste0(
      "story_local_wcb_near_threshold_shortfall()+",
      "story_local_wcb_standard_t_counterfactual()+",
      "story_local_wcb_non_near_threshold_miss_concentration()"
    )
  )
  expect_equal(ladder$threshold, 0.70, tolerance = 0)
  expect_identical(ladder$required_rejection_count, 36L)
  expect_equal(ladder$required_power_to_clear, 0.72, tolerance = 1e-12)

  repair_rungs <- ladder$repair_rungs
  expect_type(repair_rungs, "list")
  expect_identical(
    vapply(repair_rungs, `[[`, character(1L), "rung_id"),
    c(
      "observed-wild-bootstrap",
      "all-near-threshold-misses-flip",
      "standard-t-counterfactual",
      "minimum-non-near-threshold-case-combo",
      "minimum-non-near-threshold-bucket-combo"
    )
  )
  expect_identical(
    vapply(repair_rungs, `[[`, integer(1L), "rejection_count"),
    c(16L, 19L, 21L, 36L, 37L)
  )
  expect_equal(
    vapply(repair_rungs, `[[`, numeric(1L), "power"),
    c(0.32, 0.38, 0.42, 0.72, 0.74),
    tolerance = 1e-12
  )
  expect_identical(
    vapply(repair_rungs, `[[`, integer(1L), "residual_shortfall"),
    c(20L, 17L, 15L, 0L, 0L)
  )
  expect_identical(
    vapply(repair_rungs, `[[`, logical(1L), "clears_threshold"),
    c(FALSE, FALSE, FALSE, TRUE, TRUE)
  )
  expect_identical(
    repair_rungs[[4L]]$replay_seed_combo,
    c(
      58L, 72L, 60L, 91L, 84L, 69L, 85L, 45L, 83L, 82L,
      59L, 73L, 74L, 71L, 51L, 52L, 61L, 70L, 87L, 54L
    )
  )
  expect_equal(
    repair_rungs[[4L]]$share_of_non_near_threshold_misses,
    20 / 31,
    tolerance = 1e-12
  )
  expect_identical(
    repair_rungs[[5L]]$treated_cluster_counts,
    c(5L, 6L, 4L)
  )
  expect_equal(
    repair_rungs[[5L]]$share_of_non_near_threshold_misses,
    21 / 31,
    tolerance = 1e-12
  )
  expect_identical(
    ladder$summary$first_threshold_clearing_rung,
    "minimum-non-near-threshold-case-combo"
  )
  expect_identical(
    ladder$summary$minimum_threshold_clearing_case_combo,
    c(
      58L, 72L, 60L, 91L, 84L, 69L, 85L, 45L, 83L, 82L,
      59L, 73L, 74L, 71L, 51L, 52L, 61L, 70L, 87L, 54L
    )
  )
  expect_identical(
    ladder$summary$minimum_threshold_clearing_bucket_combo,
    c(5L, 6L, 4L)
  )
  expect_true(isTRUE(
    ladder$summary$minimum_case_combo_clears_with_fewer_repairs_than_bucket_combo
  ))
  expect_true(isTRUE(
    ladder$summary$near_threshold_alone_still_below_standard_t_counterfactual
  ))
  expect_true(isTRUE(
    ladder$summary$standard_t_counterfactual_still_below_threshold
  ))
  expect_true(isTRUE(
    ladder$summary$threshold_crossing_requires_non_near_threshold_combo
  ))
})

test_that("story-local TC-9.4.19 diagnostics expose non-near-threshold miss concentration", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary

  expect_type(miss_summary, "list")
  expect_identical(miss_summary$total_misses, 34L)
  expect_identical(miss_summary$near_threshold_misses, 3L)
  expect_identical(miss_summary$non_near_threshold_misses, 31L)
  expect_equal(
    miss_summary$non_near_threshold_share_of_misses,
    31 / 34,
    tolerance = 1e-12
  )
  expect_equal(
    unname(unlist(miss_summary$near_threshold_window)),
    c(0.04, 0.06),
    tolerance = 1e-12
  )

  bucket_ranking <- miss_summary$bucket_ranking
  expect_type(bucket_ranking, "list")
  expect_length(bucket_ranking, 8L)
  expect_identical(
    vapply(bucket_ranking[1:3], `[[`, integer(1L), "treated_cluster_count"),
    c(5L, 6L, 4L)
  )
  expect_identical(
    vapply(bucket_ranking[1:3], `[[`, integer(1L), "miss_count"),
    c(12L, 5L, 4L)
  )
  expect_equal(
    vapply(bucket_ranking[1:3], `[[`, numeric(1L), "avg_abs_t_stat"),
    c(1.23892960540929, 0.834258394861495, 1.22554202143167),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(bucket_ranking[1:3], `[[`, numeric(1L), "median_abs_t_stat"),
    c(1.27037497204229, 0.932053560595048, 1.18333294302385),
    tolerance = 1e-12
  )

  representative_large_misses <- lapply(
    bucket_ranking[1:3],
    `[[`,
    "representative_large_miss"
  )
  expect_identical(
    vapply(representative_large_misses, `[[`, integer(1L), "attempt_id"),
    c(24L, 15L, 36L)
  )
  expect_identical(
    vapply(representative_large_misses, `[[`, integer(1L), "replay_seed"),
    c(65L, 56L, 77L)
  )
  expect_identical(
    vapply(representative_large_misses, `[[`, integer(1L), "replay_seed_offset"),
    c(23L, 14L, 35L)
  )
  expect_equal(
    vapply(representative_large_misses, `[[`, numeric(1L), "pvalue_gap"),
    c(0.6482421875, 0.815234375, 0.643359375),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(representative_large_misses, `[[`, numeric(1L), "abs_t_stat_original"),
    c(0.420512965833072, 0.207481796923988, 0.447069219137332),
    tolerance = 1e-12
  )

  dominant_bucket <- miss_summary$dominant_bucket
  expect_type(dominant_bucket, "list")
  expect_identical(dominant_bucket$treated_cluster_count, 5L)
  expect_identical(dominant_bucket$miss_count, 12L)
  expect_equal(
    dominant_bucket$share_of_non_near_threshold_misses,
    12 / 31,
    tolerance = 1e-12
  )
  expect_equal(
    dominant_bucket$avg_pvalue,
    0.3042805989583333,
    tolerance = 1e-12
  )
  expect_equal(dominant_bucket$median_pvalue, 0.2587890625, tolerance = 0)
  expect_equal(dominant_bucket$max_pvalue, 0.6982421875, tolerance = 0)
  expect_equal(
    dominant_bucket$avg_abs_t_stat,
    1.23892960540929,
    tolerance = 1e-12
  )
  expect_equal(
    dominant_bucket$median_abs_t_stat,
    1.27037497204229,
    tolerance = 1e-12
  )
  expect_identical(dominant_bucket$representative_large_miss$attempt_id, 24L)
  expect_identical(dominant_bucket$representative_large_miss$replay_seed, 65L)
  expect_equal(
    dominant_bucket$representative_large_miss$pvalue_gap,
    0.6482421875,
    tolerance = 1e-12
  )
  expect_equal(
    dominant_bucket$representative_large_miss$abs_t_stat_original,
    0.420512965833072,
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 miss summary exposes cumulative bucket repair frontier", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary
  repair_frontier <- miss_summary$bucket_repair_frontier

  expect_type(repair_frontier, "list")
  expect_length(repair_frontier, 8L)
  expect_identical(
    vapply(repair_frontier[1:3], `[[`, integer(1L), "treated_cluster_count"),
    c(5L, 6L, 4L)
  )
  expect_identical(
    vapply(repair_frontier[1:3], `[[`, integer(1L), "bucket_rank"),
    c(1L, 2L, 3L)
  )
  expect_identical(
    vapply(repair_frontier[1:3], `[[`, integer(1L), "cumulative_miss_count"),
    c(12L, 17L, 21L)
  )
  expect_identical(
    vapply(repair_frontier[1:3], `[[`, integer(1L), "counterfactual_rejection_count"),
    c(28L, 33L, 37L)
  )
  expect_equal(
    vapply(repair_frontier[1:3], `[[`, numeric(1L), "counterfactual_power"),
    c(0.56, 0.66, 0.74),
    tolerance = 1e-12
  )
  expect_identical(
    vapply(repair_frontier[1:3], `[[`, integer(1L), "residual_shortfall"),
    c(8L, 3L, 0L)
  )
  expect_identical(
    vapply(repair_frontier[1:3], `[[`, logical(1L), "clears_threshold"),
    c(FALSE, FALSE, TRUE)
  )
  expect_equal(
    vapply(repair_frontier[1:3], `[[`, numeric(1L), "cumulative_share_of_non_near_threshold_misses"),
    c(12 / 31, 17 / 31, 21 / 31),
    tolerance = 1e-12
  )

  expect_identical(miss_summary$minimum_bucket_rank_to_clear_threshold, 3L)
  expect_identical(
    miss_summary$minimum_bucket_combo_to_clear_threshold,
    c(5L, 6L, 4L)
  )
  expect_equal(
    miss_summary$minimum_bucket_combo_share_of_non_near_threshold_misses,
    21 / 31,
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 miss summary exposes cumulative case repair frontier", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary
  case_frontier <- miss_summary$case_repair_frontier

  expect_type(case_frontier, "list")
  expect_length(case_frontier, 31L)
  expect_identical(
    vapply(case_frontier[1:5], `[[`, integer(1L), "replay_seed"),
    c(58L, 72L, 60L, 91L, 84L)
  )
  expect_identical(
    vapply(case_frontier[1:5], `[[`, integer(1L), "treated_cluster_count"),
    c(4L, 8L, 2L, 2L, 5L)
  )
  expect_equal(
    vapply(case_frontier[1:5], `[[`, numeric(1L), "pvalue_gap_above_threshold"),
    c(0.0271484375, 0.0291015625, 0.0447265625, 0.0681640625, 0.0828125),
    tolerance = 1e-15
  )
  expect_identical(
    vapply(case_frontier[1:5], `[[`, integer(1L), "counterfactual_rejection_count"),
    c(17L, 18L, 19L, 20L, 21L)
  )
  expect_identical(
    vapply(case_frontier[1:5], `[[`, integer(1L), "residual_shortfall"),
    c(19L, 18L, 17L, 16L, 15L)
  )
  expect_equal(
    vapply(case_frontier[c(1L, 5L, 20L)], `[[`, numeric(1L), "cumulative_share_of_non_near_threshold_misses"),
    c(1 / 31, 5 / 31, 20 / 31),
    tolerance = 1e-15
  )
  expect_identical(miss_summary$minimum_case_rank_to_clear_threshold, 20L)
  expect_identical(
    miss_summary$minimum_case_replay_seed_combo_to_clear_threshold,
    c(
      58L, 72L, 60L, 91L, 84L, 69L, 85L, 45L, 83L, 82L,
      59L, 73L, 74L, 71L, 51L, 52L, 61L, 70L, 87L, 54L
    )
  )
  expect_equal(
    miss_summary$minimum_case_combo_share_of_non_near_threshold_misses,
    20 / 31,
    tolerance = 1e-15
  )
})

test_that("story-local TC-9.4.19 representative large misses expose signal-vs-se ratios", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  bucket_ranking <- diagnostics$result$replication_trace$non_near_threshold_miss_summary$bucket_ranking

  top_buckets <- bucket_ranking[1:3]
  expect_identical(
    vapply(top_buckets, `[[`, integer(1L), "treated_cluster_count"),
    c(5L, 6L, 4L)
  )
  expect_equal(
    vapply(top_buckets, `[[`, numeric(1L), "avg_abs_att"),
    c(1.7035484052019, 0.743686552124342, 1.28891337279463),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(top_buckets, `[[`, numeric(1L), "avg_original_se"),
    c(1.41267080475132, 0.948591214308955, 1.07497858697201),
    tolerance = 1e-12
  )

  representative_large_misses <- lapply(
    top_buckets,
    `[[`,
    "representative_large_miss"
  )
  expect_equal(
    vapply(representative_large_misses, `[[`, numeric(1L), "abs_att"),
    c(0.697135039898459, 0.243292132423837, 0.371384544910709),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(representative_large_misses, `[[`, numeric(1L), "original_se"),
    c(1.65782055855846, 1.17259507113758, 0.830709270540556),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(representative_large_misses, `[[`, numeric(1L), "abs_att_to_bucket_avg_ratio"),
    c(0.409225260503142, 0.327143380136259, 0.288137707893798),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(representative_large_misses, `[[`, numeric(1L), "original_se_to_bucket_avg_ratio"),
    c(1.17353636316587, 1.23614371865315, 0.772768202649029),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(representative_large_misses, `[[`, numeric(1L), "abs_t_stat_to_bucket_avg_ratio"),
    c(0.339416351015482, 0.248702078638879, 0.364793055904413),
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 representative far misses classify mechanism groups", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()
  miss_summary <- diagnostics$result$replication_trace$non_near_threshold_miss_summary
  mechanism_summary <- miss_summary$representative_far_miss_mechanism_summary

  expect_type(mechanism_summary, "list")
  expect_identical(
    mechanism_summary$top_bucket_order,
    c(5L, 6L, 4L)
  )
  expect_identical(
    mechanism_summary$representative_replay_seeds,
    c(65L, 56L, 77L)
  )
  expect_identical(
    mechanism_summary$low_signal_plus_elevated_se_buckets,
    c(5L, 6L)
  )
  expect_identical(
    mechanism_summary$low_signal_dominant_buckets,
    4L
  )
  expect_identical(
    mechanism_summary$mechanism_split,
    "buckets-5-6-low-signal-plus-elevated-se_bucket-4-low-signal-dominant"
  )
  expect_true(isTRUE(mechanism_summary$all_representative_abs_att_ratios_below_one))
  expect_true(isTRUE(mechanism_summary$all_buckets_ranked_by_non_near_threshold_mass))

  bucket_mechanisms <- mechanism_summary$bucket_mechanisms
  expect_type(bucket_mechanisms, "list")
  expect_length(bucket_mechanisms, 3L)
  expect_identical(
    vapply(bucket_mechanisms, `[[`, integer(1L), "treated_cluster_count"),
    c(5L, 6L, 4L)
  )
  expect_identical(
    vapply(bucket_mechanisms, `[[`, character(1L), "mechanism_label"),
    c(
      "low-signal-plus-elevated-se",
      "low-signal-plus-elevated-se",
      "low-signal-dominant"
    )
  )
  expect_equal(
    vapply(bucket_mechanisms, `[[`, numeric(1L), "representative_abs_att_to_bucket_avg_ratio"),
    c(0.409225260503142, 0.327143380136259, 0.288137707893798),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(bucket_mechanisms, `[[`, numeric(1L), "representative_original_se_to_bucket_avg_ratio"),
    c(1.17353636316587, 1.23614371865315, 0.772768202649029),
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 diagnostics can replay the frozen near-threshold seeds", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(
    replay_seeds = c(88L, 42L, 80L, 81L, 75L)
  )
  replay_cases <- diagnostics$result$targeted_replays

  expect_type(replay_cases, "list")
  expect_length(replay_cases, 5L)
  expect_identical(
    vapply(replay_cases, `[[`, integer(1L), "replication_id"),
    c(47L, 1L, 39L, 40L, 34L)
  )
  expect_identical(
    vapply(replay_cases, `[[`, integer(1L), "seed"),
    c(88L, 42L, 80L, 81L, 75L)
  )
  expect_identical(
    vapply(replay_cases, `[[`, integer(1L), "treated_cluster_count"),
    c(6L, 8L, 7L, 5L, 4L)
  )
  expect_identical(
    vapply(replay_cases, `[[`, logical(1L), "reject"),
    c(TRUE, FALSE, FALSE, FALSE, TRUE)
  )
  expect_equal(
    vapply(replay_cases, `[[`, numeric(1L), "pvalue"),
    c(0.0478515625, 0.0576171875, 0.05859375, 0.0595703125, 0.0400390625),
    tolerance = 0
  )
  expect_equal(
    vapply(replay_cases, `[[`, numeric(1L), "att"),
    c(
      3.6154459052930346,
      4.935644839574082,
      -2.5793877511070122,
      2.814960988218609,
      2.805461571884962
    ),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(replay_cases, `[[`, numeric(1L), "original_se"),
    c(
      1.4436341616445825,
      1.2391170010075931,
      1.1323670885326091,
      1.2720609961679905,
      1.1406712782087591
    ),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(replay_cases, `[[`, numeric(1L), "t_stat_original"),
    c(
      2.5044058954481465,
      3.9831951587789063,
      -2.277872411895635,
      2.2129135290670137,
      2.4594829601482453
    ),
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 targeted replay does not depend on aggregate n_simulations", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(
    n_simulations = 1L,
    replay_seeds = 88L
  )
  replay_case <- diagnostics$result$targeted_replays[[1L]]

  expect_identical(replay_case$seed, 88L)
  expect_identical(replay_case$attempt_id, 47L)
  expect_identical(replay_case$replication_id, 47L)
  expect_identical(replay_case$treated_cluster_count, 6L)
  expect_true(isTRUE(replay_case$reject))
  expect_equal(replay_case$pvalue, 0.0478515625, tolerance = 0)
  expect_equal(replay_case$att, 3.6154459052930346, tolerance = 1e-12)
  expect_equal(replay_case$original_se, 1.4436341616445825, tolerance = 1e-12)
  expect_equal(replay_case$t_stat_original, 2.5044058954481465, tolerance = 1e-12)
})

test_that("story-local TC-9.4.19 targeted replay exposes a machine-readable consumer summary", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(
    replay_seeds = c(88L, 42L, 80L, 81L, 75L)
  )
  replay_summary <- diagnostics$result$targeted_replay_summary

  expect_type(replay_summary, "list")
  expect_identical(replay_summary$targeted_replay_count, 5L)
  expect_identical(replay_summary$replay_seed_labels, c(88L, 42L, 80L, 81L, 75L))
  expect_identical(replay_summary$replay_attempt_ids, c(47L, 1L, 39L, 40L, 34L))
  expect_identical(
    replay_summary$helper_rejection_pattern,
    c(TRUE, FALSE, FALSE, FALSE, TRUE)
  )
  expect_identical(replay_summary$treated_cluster_pattern, c(6L, 8L, 7L, 5L, 4L))
  expect_equal(
    replay_summary$weakest_hit$seed,
    88L,
    tolerance = 0
  )
  expect_equal(
    replay_summary$weakest_hit$attempt_id,
    47L,
    tolerance = 0
  )
  expect_equal(
    replay_summary$weakest_hit$pvalue,
    0.0478515625,
    tolerance = 0
  )
  expect_equal(
    replay_summary$weakest_miss$seed,
    42L,
    tolerance = 0
  )
  expect_equal(
    replay_summary$weakest_miss$attempt_id,
    1L,
    tolerance = 0
  )
  expect_equal(
    replay_summary$weakest_miss$pvalue,
    0.0576171875,
    tolerance = 0
  )
})

test_that("story-local TC-9.4.19 targeted replay summary exposes replay-seed offsets", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(
    seed = 100L,
    replay_seeds = c(100L, 146L)
  )
  replay_summary <- diagnostics$result$targeted_replay_summary

  expect_type(replay_summary, "list")
  expect_identical(replay_summary$replay_seed_labels, c(100L, 146L))
  expect_identical(replay_summary$replay_seed_offsets, c(0L, 46L))
  expect_identical(replay_summary$replay_attempt_ids, c(1L, 47L))
})

test_that("story-local TC-9.4.19 targeted replay preserves duplicate replay-seed labels", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(
    replay_seeds = c(88L, 88L, 42L)
  )
  targeted_replays <- diagnostics$result$targeted_replays
  replay_summary <- diagnostics$result$targeted_replay_summary

  expect_length(targeted_replays, 3L)
  expect_identical(vapply(targeted_replays, `[[`, integer(1L), "seed"), c(88L, 88L, 42L))
  expect_identical(vapply(targeted_replays, `[[`, integer(1L), "attempt_id"), c(47L, 47L, 1L))
  expect_identical(
    vapply(targeted_replays, `[[`, integer(1L), "replay_seed_offset"),
    c(46L, 46L, 0L)
  )
  expect_identical(vapply(targeted_replays, `[[`, logical(1L), "reject"), c(TRUE, TRUE, FALSE))
  expect_identical(replay_summary$targeted_replay_count, 3L)
  expect_identical(replay_summary$replay_seed_labels, c(88L, 88L, 42L))
  expect_identical(replay_summary$replay_seed_offsets, c(46L, 46L, 0L))
  expect_identical(replay_summary$replay_attempt_ids, c(47L, 47L, 1L))
  expect_identical(replay_summary$helper_rejection_pattern, c(TRUE, TRUE, FALSE))
})

test_that("story-local TC-9.4.19 targeted replay keeps empty consumer payload machine-readable", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(replay_seeds = integer(0))
  targeted_replays <- diagnostics$result$targeted_replays
  replay_summary <- diagnostics$result$targeted_replay_summary

  expect_identical(targeted_replays, list())
  expect_type(replay_summary, "list")
  expect_identical(replay_summary$targeted_replay_count, 0L)
  expect_identical(replay_summary$replay_seed_labels, integer(0))
  expect_identical(replay_summary$replay_seed_offsets, integer(0))
  expect_identical(replay_summary$replay_attempt_ids, integer(0))
  expect_identical(replay_summary$helper_rejection_pattern, logical(0))
  expect_identical(replay_summary$treated_cluster_pattern, integer(0))
  expect_null(replay_summary$weakest_hit)
  expect_null(replay_summary$weakest_miss)
})

test_that("story-local TC-9.4.19 targeted replay selector ranks the nearest-threshold cases", {
  selector <- lwdid:::.story_local_wcb_targeted_replay_selector()

  expect_identical(selector$case_id, "TC-9.4.19")
  expect_identical(
    selector$exact_status,
    "story-local-targeted-replay-selector-contract-frozen"
  )
  expect_identical(
    selector$numeric_status,
    "targeted-replay-attempt-mapping-machine-readable"
  )
  expect_identical(
    selector$blocker_boundary,
    "story-local-hardening-power-gap-tc-9-4-19"
  )
  expect_identical(
    selector$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(selector$base_seed, 42L)
  expect_equal(selector$pvalue_threshold, 0.05, tolerance = 0)
  expect_identical(selector$requested_n_bootstrap, 199L)
  expect_identical(selector$actual_n_bootstrap, 1024L)
  expect_identical(selector$replay_seed_labels, c(88L, 42L, 80L, 81L, 75L))
  expect_identical(
    vapply(selector$targeted_replays, `[[`, integer(1L), "attempt_id"),
    c(47L, 1L, 39L, 40L, 34L)
  )
  expect_identical(
    selector$summary$replay_seed_labels,
    c(88L, 42L, 80L, 81L, 75L)
  )
  expect_identical(
    selector$summary$replay_attempt_ids,
    c(47L, 1L, 39L, 40L, 34L)
  )
  expect_identical(
    selector$summary$helper_rejection_pattern,
    c(TRUE, FALSE, FALSE, FALSE, TRUE)
  )
  expect_identical(
    selector$summary$treated_cluster_pattern,
    c(6L, 8L, 7L, 5L, 4L)
  )
  expect_identical(selector$summary$weakest_hit$seed, 88L)
  expect_identical(selector$summary$weakest_miss$seed, 42L)
})

test_that("story-local TC-9.4.19 standard-t counterfactual helper freezes the paired rescue gap", {
  counterfactual <- lwdid:::.story_local_wcb_standard_t_counterfactual()

  expect_identical(counterfactual$case_id, "TC-9.4.19")
  expect_identical(
    counterfactual$exact_status,
    "story-local-standard-t-counterfactual-frozen"
  )
  expect_identical(
    counterfactual$numeric_status,
    "paired-standard-t-vs-wcb-rescue-set-machine-readable"
  )
  expect_identical(
    counterfactual$consumer_summary_source,
    "standard_t_counterfactual_summary"
  )
  expect_identical(counterfactual$paired_replication_count, 50L)
  expect_identical(counterfactual$standard_t_rejection_count, 21L)
  expect_identical(counterfactual$wild_bootstrap_rejection_count, 16L)
  expect_identical(counterfactual$shared_rejection_count, 16L)
  expect_identical(counterfactual$standard_only_rejection_count, 5L)
  expect_identical(counterfactual$wcb_only_rejection_count, 0L)
  expect_identical(
    counterfactual$standard_only_replay_seeds,
    c(42L, 60L, 74L, 80L, 91L)
  )
  expect_identical(
    vapply(
      counterfactual$standard_only_bucket_ranking,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ),
    c(2L, 7L, 8L, 9L)
  )
  expect_identical(
    vapply(
      counterfactual$standard_only_repair_frontier,
      `[[`,
      integer(1L),
      "replay_seed"
    ),
    c(42L, 80L, 60L, 91L, 74L)
  )
  expect_true(isTRUE(
    counterfactual$summary_consistency_checks$standard_only_seed_vector_matches_cases
  ))
  expect_true(isTRUE(
    counterfactual$summary_consistency_checks$standard_only_bucket_ranking_matches_cases
  ))
  expect_true(isTRUE(
    counterfactual$summary_consistency_checks$standard_only_cases_match_replications
  ))
  expect_true(isTRUE(
    counterfactual$summary_consistency_checks$wcb_only_cases_remain_empty
  ))
  expect_equal(counterfactual$summary$standard_t_power, 0.42, tolerance = 1e-12)
  expect_equal(
    counterfactual$summary$wild_bootstrap_power,
    0.32,
    tolerance = 1e-12
  )
  expect_equal(
    counterfactual$summary$rejection_gap_share_of_replications,
    0.1,
    tolerance = 1e-12
  )
  expect_identical(counterfactual$summary$closest_standard_only_seed_to_threshold, 42L)
  expect_identical(
    counterfactual$summary$closest_non_near_standard_only_seed,
    60L
  )
})

test_that("story-local TC-9.4.19 standard-t rescue threshold split helper partitions near-threshold rescues", {
  threshold_split <- lwdid:::.story_local_wcb_standard_t_rescue_threshold_split()

  expect_identical(threshold_split$case_id, "TC-9.4.19")
  expect_identical(
    threshold_split$exact_status,
    "story-local-standard-t-rescue-threshold-split-frozen"
  )
  expect_identical(
    threshold_split$numeric_status,
    "standard-only-rescues-threshold-split-machine-readable"
  )
  expect_identical(
    threshold_split$consumer_summary_source,
    "standard_t_counterfactual_summary"
  )
  expect_identical(threshold_split$standard_only_rejection_count, 5L)
  expect_identical(threshold_split$standard_only_near_threshold_count, 2L)
  expect_identical(threshold_split$standard_only_non_near_threshold_count, 3L)
  expect_identical(
    threshold_split$standard_only_near_threshold_replay_seeds,
    c(42L, 80L)
  )
  expect_identical(
    threshold_split$standard_only_non_near_threshold_replay_seeds,
    c(60L, 74L, 91L)
  )
  expect_identical(
    vapply(
      threshold_split$standard_only_near_threshold_cases,
      `[[`,
      integer(1L),
      "replay_seed"
    ),
    c(42L, 80L)
  )
  expect_identical(
    vapply(
      threshold_split$standard_only_non_near_threshold_cases,
      `[[`,
      integer(1L),
      "replay_seed"
    ),
    c(60L, 74L, 91L)
  )
  expect_true(isTRUE(
    threshold_split$summary_consistency_checks$standard_only_partition_is_exhaustive
  ))
  expect_true(isTRUE(
    threshold_split$summary_consistency_checks$near_threshold_cases_match_window
  ))
  expect_true(isTRUE(
    threshold_split$summary_consistency_checks$non_near_threshold_cases_exceed_window
  ))
  expect_equal(
    threshold_split$summary$near_threshold_share_of_standard_only,
    0.4,
    tolerance = 1e-12
  )
  expect_equal(
    threshold_split$summary$non_near_threshold_share_of_standard_only,
    0.6,
    tolerance = 1e-12
  )
  expect_identical(
    threshold_split$summary$closest_standard_only_seed_to_threshold,
    42L
  )
  expect_identical(
    threshold_split$summary$closest_non_near_standard_only_seed,
    60L
  )
  expect_equal(
    threshold_split$standard_only_near_threshold_cases[[1L]]$wcb_gap_above_threshold,
    0.0076171875,
    tolerance = 1e-15
  )
  expect_equal(
    threshold_split$standard_only_non_near_threshold_cases[[1L]]$wcb_gap_above_window_upper,
    0.0347265625,
    tolerance = 1e-15
  )
})

test_that("story-local TC-9.4.19 near-threshold shortfall helper caps threshold-only rescue capacity", {
  shortfall <- lwdid:::.story_local_wcb_near_threshold_shortfall()

  expect_identical(shortfall$case_id, "TC-9.4.19")
  expect_identical(
    shortfall$exact_status,
    "story-local-near-threshold-shortfall-bridge-frozen"
  )
  expect_identical(
    shortfall$numeric_status,
    "near-threshold-flip-capacity-machine-readable"
  )
  expect_equal(shortfall$pvalue_threshold, 0.05, tolerance = 0)
  expect_equal(shortfall$near_threshold_window, c(0.04, 0.06), tolerance = 1e-12)
  expect_identical(shortfall$actual_rejection_count, 16L)
  expect_identical(shortfall$required_rejection_count, 36L)
  expect_identical(shortfall$rejection_count_shortfall, 20L)
  expect_identical(shortfall$near_threshold_cases$total, 5L)
  expect_identical(shortfall$near_threshold_cases$hits, 2L)
  expect_identical(shortfall$near_threshold_cases$misses, 3L)
  expect_identical(shortfall$counterfactual_if_all_near_threshold_misses_flip, 19L)
  expect_equal(
    shortfall$counterfactual_power_if_all_near_threshold_misses_flip,
    0.38,
    tolerance = 1e-12
  )
  expect_identical(
    shortfall$residual_shortfall_after_flipping_near_threshold_misses,
    17L
  )
  expect_false(isTRUE(shortfall$can_clear_threshold_via_near_threshold_flips_alone))
})

test_that("story-local TC-9.4.19 standard-t repair frontier helper replays the ordered rescue seeds", {
  frontier_replays <- lwdid:::.story_local_wcb_standard_t_repair_frontier_replays()

  expect_identical(frontier_replays$case_id, "TC-9.4.19")
  expect_identical(
    frontier_replays$exact_status,
    "story-local-standard-t-repair-frontier-replays-frozen"
  )
  expect_identical(
    frontier_replays$numeric_status,
    "repair-frontier-targeted-replays-machine-readable"
  )
  expect_identical(
    frontier_replays$consumer_summary_source,
    "standard_t_counterfactual_summary"
  )
  expect_identical(
    frontier_replays$replay_selector_source,
    "story_local_wcb_power_diagnostics(replay_seeds = standard_only_repair_frontier)"
  )
  expect_identical(
    vapply(
      frontier_replays$standard_only_repair_frontier,
      `[[`,
      integer(1L),
      "replay_seed"
    ),
    c(42L, 80L, 60L, 91L, 74L)
  )
  expect_identical(
    vapply(frontier_replays$targeted_replays, `[[`, integer(1L), "seed"),
    c(42L, 80L, 60L, 91L, 74L)
  )
  expect_identical(
    vapply(frontier_replays$targeted_replays, `[[`, logical(1L), "reject"),
    c(FALSE, FALSE, FALSE, FALSE, FALSE)
  )
  expect_true(isTRUE(
    frontier_replays$summary_consistency_checks$replay_seed_order_matches_frontier
  ))
  expect_true(isTRUE(
    frontier_replays$summary_consistency_checks$all_targeted_replays_remain_non_rejections
  ))
  expect_equal(
    unlist(frontier_replays$summary$frontier_wcb_gaps, use.names = FALSE),
    c(0.0076171875, 0.00859375, 0.0447265625, 0.0681640625, 0.2166015625),
    tolerance = 1e-15
  )
  expect_identical(frontier_replays$summary$closest_frontier_seed, 42L)
  expect_identical(frontier_replays$summary$furthest_frontier_seed, 74L)
})

test_that("story-local TC-9.4.19 non-near-threshold miss helper freezes dominant miss mass", {
  miss_concentration <- lwdid:::.story_local_wcb_non_near_threshold_miss_concentration()

  expect_identical(miss_concentration$case_id, "TC-9.4.19")
  expect_identical(
    miss_concentration$exact_status,
    "story-local-non-near-threshold-miss-concentration-frozen"
  )
  expect_identical(
    miss_concentration$numeric_status,
    "non-near-threshold-miss-mass-machine-readable"
  )
  expect_identical(miss_concentration$total_misses, 34L)
  expect_identical(miss_concentration$near_threshold_misses, 3L)
  expect_identical(miss_concentration$non_near_threshold_misses, 31L)
  expect_equal(
    miss_concentration$non_near_threshold_share_of_misses,
    31 / 34,
    tolerance = 1e-12
  )
  expect_equal(
    unname(unlist(miss_concentration$near_threshold_window)),
    c(0.04, 0.06),
    tolerance = 1e-12
  )
  expect_identical(
    vapply(
      miss_concentration$bucket_ranking[1:3],
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ),
    c(5L, 6L, 4L)
  )
  expect_identical(
    vapply(
      miss_concentration$bucket_ranking[1:3],
      `[[`,
      integer(1L),
      "miss_count"
    ),
    c(12L, 5L, 4L)
  )
  expect_identical(miss_concentration$dominant_bucket$treated_cluster_count, 5L)
  expect_identical(miss_concentration$dominant_bucket$miss_count, 12L)
  expect_equal(
    miss_concentration$dominant_bucket$avg_pvalue,
    0.3042805989583333,
    tolerance = 1e-12
  )
})

test_that("story-local TC-9.4.19 non-near-threshold case frontier helper freezes the minimum case combo", {
  case_frontier <- lwdid:::.story_local_wcb_non_near_threshold_case_frontier()

  expect_identical(case_frontier$case_id, "TC-9.4.19")
  expect_identical(
    case_frontier$exact_status,
    "story-local-non-near-threshold-case-frontier-frozen"
  )
  expect_identical(
    case_frontier$numeric_status,
    "case-level-threshold-clearing-frontier-machine-readable"
  )
  expect_identical(
    case_frontier$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(
    case_frontier$consumer_summary_source,
    "non_near_threshold_miss_summary"
  )
  expect_identical(case_frontier$non_near_threshold_misses, 31L)
  expect_identical(case_frontier$minimum_case_rank_to_clear_threshold, 20L)
  expect_identical(
    case_frontier$minimum_case_replay_seed_combo_to_clear_threshold,
    c(
      58L, 72L, 60L, 91L, 84L, 69L, 85L, 45L, 83L, 82L,
      59L, 73L, 74L, 71L, 51L, 52L, 61L, 70L, 87L, 54L
    )
  )
  expect_equal(
    case_frontier$minimum_case_combo_share_of_non_near_threshold_misses,
    20 / 31,
    tolerance = 1e-15
  )
  expect_identical(
    case_frontier$summary$minimum_case_combo_bucket_support,
    c(2L, 3L, 4L, 5L, 6L, 8L, 9L)
  )
  expect_identical(
    names(case_frontier$summary$minimum_case_combo_bucket_frequency),
    c("5", "6", "2", "4", "3", "8", "9")
  )
  expect_identical(
    unname(case_frontier$summary$minimum_case_combo_bucket_frequency),
    c(10L, 3L, 2L, 2L, 1L, 1L, 1L)
  )
  expect_true(isTRUE(
    case_frontier$summary$minimum_case_combo_spans_more_buckets_than_bucket_combo
  ))
  expect_true(isTRUE(
    case_frontier$summary$minimum_case_combo_uses_bucket_outside_bucket_level_combo
  ))
  expect_identical(case_frontier$summary$furthest_case_seed_in_minimum_combo, 54L)
  expect_equal(
    case_frontier$summary$furthest_case_gap_above_threshold_in_minimum_combo,
    0.354296875,
    tolerance = 1e-15
  )
})

test_that("story-local TC-9.4.19 case-frontier tail helper shows the crossing edge is still all misses", {
  tail_replays <- lwdid:::.story_local_wcb_non_near_threshold_case_frontier_tail_replays()

  expect_identical(tail_replays$case_id, "TC-9.4.19")
  expect_identical(
    tail_replays$exact_status,
    "story-local-non-near-threshold-case-frontier-tail-replays-frozen"
  )
  expect_identical(
    tail_replays$numeric_status,
    "threshold-edge-case-tail-replays-machine-readable"
  )
  expect_identical(
    tail_replays$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(
    tail_replays$consumer_summary_source,
    "story_local_wcb_non_near_threshold_case_frontier()"
  )
  expect_identical(
    tail_replays$replay_selector_source,
    "story_local_wcb_power_diagnostics(replay_seeds = threshold_edge_case_tail)"
  )
  expect_identical(
    tail_replays$tail_case_replay_seeds,
    c(52L, 61L, 70L, 87L, 54L)
  )

  tail_cases <- tail_replays$tail_case_frontier
  expect_type(tail_cases, "list")
  expect_length(tail_cases, 5L)
  expect_identical(
    vapply(tail_cases, `[[`, integer(1L), "case_rank"),
    16:20
  )
  expect_identical(
    vapply(tail_cases, `[[`, integer(1L), "replay_seed"),
    c(52L, 61L, 70L, 87L, 54L)
  )
  expect_identical(
    vapply(tail_cases, `[[`, integer(1L), "treated_cluster_count"),
    c(6L, 5L, 5L, 3L, 6L)
  )
  expect_equal(
    vapply(tail_cases, `[[`, numeric(1L), "pvalue_gap_above_threshold"),
    c(0.2859375, 0.3015625, 0.3337890625, 0.34453125, 0.354296875),
    tolerance = 1e-15
  )
  expect_identical(
    vapply(tail_cases, `[[`, integer(1L), "counterfactual_rejection_count"),
    c(32L, 33L, 34L, 35L, 36L)
  )
  expect_equal(
    vapply(tail_cases, `[[`, numeric(1L), "counterfactual_power"),
    c(0.64, 0.66, 0.68, 0.70, 0.72),
    tolerance = 1e-12
  )
  expect_identical(
    vapply(tail_cases, `[[`, integer(1L), "residual_shortfall"),
    c(4L, 3L, 2L, 1L, 0L)
  )
  expect_identical(
    vapply(
      lapply(tail_cases, `[[`, "targeted_replay"),
      `[[`,
      integer(1L),
      "seed"
    ),
    c(52L, 61L, 70L, 87L, 54L)
  )
  expect_identical(
    vapply(
      lapply(tail_cases, `[[`, "targeted_replay"),
      `[[`,
      logical(1L),
      "reject"
    ),
    c(FALSE, FALSE, FALSE, FALSE, FALSE)
  )

  expect_true(isTRUE(tail_replays$summary$all_tail_targeted_replays_remain_misses))
  expect_identical(tail_replays$summary$pre_threshold_case_rank, 19L)
  expect_identical(tail_replays$summary$pre_threshold_replay_seed, 87L)
  expect_equal(
    tail_replays$summary$pre_threshold_counterfactual_power,
    0.70,
    tolerance = 1e-12
  )
  expect_identical(tail_replays$summary$threshold_crossing_case_rank, 20L)
  expect_identical(tail_replays$summary$threshold_crossing_replay_seed, 54L)
  expect_true(isTRUE(tail_replays$summary$threshold_crossing_seed_still_observed_miss))
  expect_identical(tail_replays$summary$lowest_gap_tail_seed, 52L)
  expect_true(isTRUE(tail_replays$summary$threshold_crossing_seed_has_largest_gap_in_tail))
})

test_that("story-local TC-9.4.19 residual case helper isolates the harder post-threshold miss family", {
  residual_replays <-
    lwdid:::.story_local_wcb_non_near_threshold_case_frontier_residual_replays()

  expect_identical(residual_replays$case_id, "TC-9.4.19")
  expect_identical(
    residual_replays$exact_status,
    "story-local-non-near-threshold-case-frontier-residual-replays-frozen"
  )
  expect_identical(
    residual_replays$numeric_status,
    "post-threshold-residual-replays-machine-readable"
  )
  expect_identical(
    residual_replays$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(
    residual_replays$consumer_summary_source,
    "story_local_wcb_non_near_threshold_case_frontier()"
  )
  expect_identical(
    residual_replays$replay_selector_source,
    "story_local_wcb_power_diagnostics(replay_seeds = residual_case_frontier)"
  )
  expect_identical(
    residual_replays$residual_case_replay_seeds,
    c(57L, 48L, 44L, 86L, 55L, 77L, 65L, 64L, 56L, 53L, 79L)
  )

  residual_cases <- residual_replays$residual_case_frontier
  expect_type(residual_cases, "list")
  expect_length(residual_cases, 11L)
  expect_identical(
    vapply(residual_cases, `[[`, integer(1L), "case_rank"),
    21:31
  )
  expect_identical(
    vapply(residual_cases, `[[`, integer(1L), "replay_seed"),
    c(57L, 48L, 44L, 86L, 55L, 77L, 65L, 64L, 56L, 53L, 79L)
  )
  expect_identical(
    vapply(residual_cases, `[[`, integer(1L), "treated_cluster_count"),
    c(2L, 7L, 3L, 4L, 5L, 4L, 5L, 6L, 6L, 7L, 7L)
  )
  expect_equal(
    vapply(residual_cases, `[[`, numeric(1L), "pvalue_gap_above_threshold"),
    c(
      0.373828125,
      0.43828125,
      0.440234375,
      0.528125,
      0.5330078125,
      0.643359375,
      0.6482421875,
      0.6990234375,
      0.815234375,
      0.8435546875,
      0.89921875
    ),
    tolerance = 1e-15
  )
  expect_identical(
    vapply(
      lapply(residual_cases, `[[`, "targeted_replay"),
      `[[`,
      integer(1L),
      "seed"
    ),
    c(57L, 48L, 44L, 86L, 55L, 77L, 65L, 64L, 56L, 53L, 79L)
  )
  expect_identical(
    vapply(
      lapply(residual_cases, `[[`, "targeted_replay"),
      `[[`,
      logical(1L),
      "reject"
    ),
    rep(FALSE, 11L)
  )

  expect_identical(residual_replays$summary$residual_case_count, 11L)
  expect_identical(residual_replays$summary$threshold_clearing_rank, 20L)
  expect_identical(residual_replays$summary$first_residual_case_rank, 21L)
  expect_identical(residual_replays$summary$easiest_residual_replay_seed, 57L)
  expect_equal(
    residual_replays$summary$easiest_residual_gap_above_threshold,
    0.373828125,
    tolerance = 1e-15
  )
  expect_identical(residual_replays$summary$hardest_residual_replay_seed, 79L)
  expect_equal(
    residual_replays$summary$hardest_residual_gap_above_threshold,
    0.89921875,
    tolerance = 1e-15
  )
  expect_identical(
    residual_replays$summary$residual_bucket_support,
    c(2L, 3L, 4L, 5L, 6L, 7L)
  )
  expect_identical(
    names(residual_replays$summary$residual_bucket_frequency),
    c("7", "4", "5", "6", "2", "3")
  )
  expect_identical(
    unname(residual_replays$summary$residual_bucket_frequency),
    c(3L, 2L, 2L, 2L, 1L, 1L)
  )
  expect_true(isTRUE(
    residual_replays$summary$all_residual_targeted_replays_remain_misses
  ))
  expect_true(isTRUE(
    residual_replays$summary$bucket7_only_enters_after_threshold
  ))
  expect_true(isTRUE(
    residual_replays$summary$residual_family_excludes_bucket8_and_9
  ))
})

test_that("story-local TC-9.4.19 threshold-crossing partition keeps the residual search space machine-readable", {
  partition <- lwdid:::.story_local_wcb_case_combo_residual_partition()

  expect_identical(partition$case_id, "TC-9.4.19")
  expect_identical(
    partition$exact_status,
    "story-local-case-combo-residual-partition-frozen"
  )
  expect_identical(
    partition$numeric_status,
    "threshold-clearing-vs-residual-frontier-machine-readable"
  )
  expect_identical(
    partition$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(
    partition$consumer_summary_source,
    paste0(
      "story_local_wcb_non_near_threshold_case_frontier()+",
      "story_local_wcb_non_near_threshold_case_frontier_residual_replays()"
    )
  )
  expect_identical(partition$minimum_case_rank_to_clear_threshold, 20L)
  expect_identical(
    partition$minimum_case_combo_replay_seeds,
    c(
      58L, 72L, 60L, 91L, 84L, 69L, 85L, 45L, 83L, 82L,
      59L, 73L, 74L, 71L, 51L, 52L, 61L, 70L, 87L, 54L
    )
  )
  expect_identical(
    partition$residual_case_replay_seeds,
    c(57L, 48L, 44L, 86L, 55L, 77L, 65L, 64L, 56L, 53L, 79L)
  )
  expect_identical(partition$overlap_bucket_support, c(2L, 3L, 4L, 5L, 6L))
  expect_identical(partition$minimum_only_bucket_support, c(8L, 9L))
  expect_identical(partition$residual_only_bucket_support, 7L)
  expect_identical(partition$summary$non_near_threshold_miss_count, 31L)
  expect_identical(partition$summary$minimum_case_combo_count, 20L)
  expect_identical(partition$summary$residual_case_count, 11L)
  expect_equal(
    partition$summary$minimum_case_combo_share_of_non_near_threshold_misses,
    20 / 31,
    tolerance = 1e-15
  )
  expect_equal(
    partition$summary$residual_case_share_of_non_near_threshold_misses,
    11 / 31,
    tolerance = 1e-15
  )
  expect_true(isTRUE(
    partition$summary$combined_partition_reconstructs_non_near_threshold_misses
  ))
  expect_true(isTRUE(
    partition$summary$residual_frontier_starts_after_threshold
  ))
  expect_true(isTRUE(partition$summary$bucket7_is_residual_only))
  expect_true(isTRUE(partition$summary$buckets8_and_9_are_threshold_only))
})

test_that("story-local TC-9.4.19 top-bucket replay helper freezes representative far misses", {
  top_bucket_replays <- lwdid:::.story_local_wcb_top_bucket_representative_replays()

  expect_identical(top_bucket_replays$case_id, "TC-9.4.19")
  expect_identical(
    top_bucket_replays$exact_status,
    "story-local-top-bucket-representative-replays-frozen"
  )
  expect_identical(
    top_bucket_replays$numeric_status,
    "top-bucket-far-miss-replays-machine-readable"
  )
  expect_identical(top_bucket_replays$total_misses, 34L)
  expect_identical(top_bucket_replays$non_near_threshold_misses, 31L)
  expect_identical(top_bucket_replays$top_bucket_order, c(5L, 6L, 4L))
  expect_identical(top_bucket_replays$representative_replay_seeds, c(65L, 56L, 77L))
  expect_equal(
    top_bucket_replays$summary$top_three_share_of_non_near_threshold_misses,
    21 / 31,
    tolerance = 1e-12
  )
  expect_equal(
    top_bucket_replays$summary$bucket5_share_of_top_three_misses,
    12 / 21,
    tolerance = 1e-12
  )
  expect_identical(
    top_bucket_replays$summary$bucket5_excess_miss_count_vs_bucket6,
    7L
  )
  expect_identical(
    top_bucket_replays$summary$bucket5_excess_miss_count_vs_bucket4,
    8L
  )
  expect_true(top_bucket_replays$summary$all_replayed_cases_remain_misses)
  expect_identical(
    vapply(
      top_bucket_replays$top_bucket_replays,
      `[[`,
      integer(1L),
      "treated_cluster_count"
    ),
    c(5L, 6L, 4L)
  )
  expect_identical(
    vapply(
      top_bucket_replays$top_bucket_replays,
      `[[`,
      integer(1L),
      "miss_count"
    ),
    c(12L, 5L, 4L)
  )
  expect_identical(
    vapply(
      lapply(
        top_bucket_replays$top_bucket_replays,
        `[[`,
        "representative_large_miss"
      ),
      `[[`,
      integer(1L),
      "replay_seed"
    ),
    c(65L, 56L, 77L)
  )
  expect_equal(
    vapply(
      lapply(top_bucket_replays$top_bucket_replays, `[[`, "targeted_replay"),
      `[[`,
      numeric(1L),
      "pvalue"
    ),
    c(0.6982421875, 0.865234375, 0.693359375),
    tolerance = 1e-15
  )
  expect_true(all(!vapply(
    lapply(top_bucket_replays$top_bucket_replays, `[[`, "targeted_replay"),
    `[[`,
    logical(1L),
    "reject"
  )))
})

test_that("story-local TC-9.4.19 representative far-miss mechanism helper exposes the bucket split", {
  mechanism_summary <- lwdid:::.story_local_wcb_representative_far_miss_mechanism_summary()

  expect_identical(mechanism_summary$case_id, "TC-9.4.19")
  expect_identical(
    mechanism_summary$exact_status,
    "story-local-representative-far-miss-mechanism-frozen"
  )
  expect_identical(
    mechanism_summary$numeric_status,
    "representative-far-miss-signal-vs-se-mechanism-machine-readable"
  )
  expect_identical(
    mechanism_summary$consumer_summary_source,
    "representative_far_miss_mechanism_summary"
  )
  expect_identical(mechanism_summary$top_bucket_order, c(5L, 6L, 4L))
  expect_identical(mechanism_summary$representative_replay_seeds, c(65L, 56L, 77L))
  expect_identical(
    mechanism_summary$low_signal_plus_elevated_se_buckets,
    c(5L, 6L)
  )
  expect_identical(mechanism_summary$low_signal_dominant_buckets, 4L)
  expect_identical(
    mechanism_summary$summary$mechanism_split,
    "buckets-5-6-low-signal-plus-elevated-se_bucket-4-low-signal-dominant"
  )
  expect_true(
    mechanism_summary$summary_consistency_checks$top_bucket_order_matches_bucket_ranking
  )
  expect_true(
    mechanism_summary$summary_consistency_checks$representative_replay_seeds_match_bucket_ranking
  )
  expect_true(
    mechanism_summary$summary_consistency_checks$mechanism_labels_match_bucket_ranking
  )
  expect_true(
    mechanism_summary$summary_consistency_checks$ratios_match_bucket_ranking
  )
  expect_identical(
    vapply(
      mechanism_summary$bucket_mechanisms,
      `[[`,
      character(1L),
      "mechanism_label"
    ),
    c(
      "low-signal-plus-elevated-se",
      "low-signal-plus-elevated-se",
      "low-signal-dominant"
    )
  )
  expect_equal(
    mechanism_summary$bucket_mechanisms[[1L]]$representative_abs_att_to_bucket_avg_ratio,
    0.409225260503,
    tolerance = 1e-11
  )
  expect_equal(
    mechanism_summary$bucket_mechanisms[[1L]]$representative_original_se_to_bucket_avg_ratio,
    1.17353636317,
    tolerance = 1e-11
  )
  expect_equal(
    mechanism_summary$bucket_mechanisms[[3L]]$representative_abs_att_to_bucket_avg_ratio,
    0.288137707894,
    tolerance = 1e-11
  )
  expect_equal(
    mechanism_summary$bucket_mechanisms[[3L]]$representative_original_se_to_bucket_avg_ratio,
    0.772768202649,
    tolerance = 1e-11
  )
})

test_that("story-local TC-9.4.19 case-combo mechanism bridge quantifies which bucket families dominate the threshold-clearing combo", {
  bridge <- lwdid:::.story_local_wcb_case_combo_mechanism_bridge()

  expect_identical(bridge$case_id, "TC-9.4.19")
  expect_identical(
    bridge$exact_status,
    "story-local-case-combo-mechanism-bridge-frozen"
  )
  expect_identical(
    bridge$numeric_status,
    "case-combo-mechanism-family-shares-machine-readable"
  )
  expect_identical(
    bridge$consumer_summary_source,
    paste0(
      "story_local_wcb_non_near_threshold_case_frontier()+",
      "story_local_wcb_representative_far_miss_mechanism_summary()"
    )
  )
  expect_identical(bridge$minimum_case_rank_to_clear_threshold, 20L)
  expect_identical(bridge$top_bucket_order, c(5L, 6L, 4L))
  expect_identical(
    bridge$low_signal_plus_elevated_se_buckets,
    c(5L, 6L)
  )
  expect_identical(bridge$low_signal_dominant_buckets, 4L)
  expect_identical(bridge$residual_bucket_support, c(2L, 3L, 8L, 9L))
  expect_identical(
    names(bridge$minimum_case_combo_bucket_frequency),
    c("5", "6", "2", "4", "3", "8", "9")
  )
  expect_identical(
    unname(bridge$minimum_case_combo_bucket_frequency),
    c(10L, 3L, 2L, 2L, 1L, 1L, 1L)
  )
  expect_identical(bridge$summary$top_bucket_case_count, 15L)
  expect_equal(bridge$summary$top_bucket_case_share, 0.75, tolerance = 1e-12)
  expect_identical(
    bridge$summary$low_signal_plus_elevated_se_case_count,
    13L
  )
  expect_equal(
    bridge$summary$low_signal_plus_elevated_se_case_share,
    0.65,
    tolerance = 1e-12
  )
  expect_identical(bridge$summary$low_signal_dominant_case_count, 2L)
  expect_equal(
    bridge$summary$low_signal_dominant_case_share,
    0.10,
    tolerance = 1e-12
  )
  expect_identical(bridge$summary$residual_case_count, 5L)
  expect_equal(bridge$summary$residual_case_share, 0.25, tolerance = 1e-12)
  expect_identical(bridge$summary$dominant_bucket_case_count, 10L)
  expect_equal(
    bridge$summary$dominant_bucket_share_of_minimum_case_combo,
    0.50,
    tolerance = 1e-12
  )
  expect_true(isTRUE(
    bridge$summary$top_bucket_family_covers_threshold_combo_majority
  ))
  expect_true(isTRUE(
    bridge$summary_consistency_checks$top_bucket_and_residual_reconstruct_case_combo
  ))
  expect_true(isTRUE(
    bridge$summary_consistency_checks$mechanism_split_reconstructs_top_bucket_case_count
  ))
})

test_that("story-local TC-9.4.19 partition-exclusive bucket helper isolates the 7-vs-8/9 replay families", {
  exclusive <- lwdid:::.story_local_wcb_partition_exclusive_bucket_representatives()

  expect_identical(exclusive$case_id, "TC-9.4.19")
  expect_identical(
    exclusive$exact_status,
    "story-local-partition-exclusive-bucket-representatives-frozen"
  )
  expect_identical(
    exclusive$numeric_status,
    "threshold-only-vs-residual-only-representative-replays-machine-readable"
  )
  expect_identical(
    exclusive$consumer_summary_source,
    paste0(
      "story_local_wcb_case_combo_residual_partition()+",
      "story_local_wcb_power_diagnostics(replay_seeds = partition_exclusive_bucket_representatives)"
    )
  )
  expect_identical(exclusive$threshold_only_bucket_support, c(8L, 9L))
  expect_identical(exclusive$residual_only_bucket_support, 7L)
  expect_identical(exclusive$threshold_only_replay_seeds, c(72L, 74L))
  expect_identical(exclusive$residual_only_replay_seeds, c(48L, 53L, 79L))

  expect_identical(
    vapply(exclusive$threshold_only_representatives, `[[`, integer(1L), "case_rank"),
    c(2L, 13L)
  )
  expect_identical(
    vapply(exclusive$threshold_only_representatives, `[[`, integer(1L), "treated_cluster_count"),
    c(8L, 9L)
  )
  expect_identical(
    vapply(exclusive$residual_only_representatives, `[[`, integer(1L), "case_rank"),
    c(22L, 30L, 31L)
  )
  expect_identical(
    vapply(exclusive$residual_only_representatives, `[[`, integer(1L), "treated_cluster_count"),
    c(7L, 7L, 7L)
  )
  expect_identical(
    vapply(
      lapply(exclusive$threshold_only_representatives, `[[`, "targeted_replay"),
      `[[`,
      integer(1L),
      "seed"
    ),
    c(72L, 74L)
  )
  expect_identical(
    vapply(
      lapply(exclusive$residual_only_representatives, `[[`, "targeted_replay"),
      `[[`,
      integer(1L),
      "seed"
    ),
    c(48L, 53L, 79L)
  )
  expect_equal(
    vapply(
      lapply(exclusive$threshold_only_representatives, `[[`, "targeted_replay"),
      `[[`,
      numeric(1L),
      "pvalue"
    ),
    c(0.0791015625, 0.2666015625),
    tolerance = 1e-15
  )
  expect_equal(
    vapply(
      lapply(exclusive$residual_only_representatives, `[[`, "targeted_replay"),
      `[[`,
      numeric(1L),
      "pvalue"
    ),
    c(0.48828125, 0.8935546875, 0.94921875),
    tolerance = 1e-15
  )

  expect_identical(exclusive$summary$threshold_only_case_count, 2L)
  expect_identical(exclusive$summary$residual_only_case_count, 3L)
  expect_equal(
    exclusive$summary$threshold_only_share_of_partition_exclusive_cases,
    0.4,
    tolerance = 1e-12
  )
  expect_equal(
    exclusive$summary$residual_only_share_of_partition_exclusive_cases,
    0.6,
    tolerance = 1e-12
  )
  expect_identical(exclusive$summary$closest_threshold_only_replay_seed, 72L)
  expect_equal(
    exclusive$summary$closest_threshold_only_gap_above_threshold,
    0.0291015625,
    tolerance = 1e-15
  )
  expect_identical(exclusive$summary$highest_threshold_only_t_stat_seed, 74L)
  expect_equal(
    exclusive$summary$highest_threshold_only_t_stat,
    6.464390434806528,
    tolerance = 1e-12
  )
  expect_identical(exclusive$summary$easiest_residual_only_replay_seed, 48L)
  expect_identical(exclusive$summary$hardest_residual_only_replay_seed, 79L)
  expect_true(isTRUE(exclusive$summary$all_partition_exclusive_targeted_replays_remain_misses))
  expect_true(isTRUE(exclusive$summary$bucket8_is_closest_threshold_only_case))
  expect_true(isTRUE(exclusive$summary$bucket9_has_highest_threshold_only_t_stat))
  expect_true(isTRUE(exclusive$summary$bucket7_reconstructs_residual_only_family))
})

test_that("story-local TC-9.4.19 partition-exclusive signal contrast keeps threshold-only and residual-only families machine-readable", {
  contrast <- lwdid:::.story_local_wcb_partition_exclusive_bucket_signal_contrast()

  expect_identical(contrast$case_id, "TC-9.4.19")
  expect_identical(
    contrast$exact_status,
    "story-local-partition-exclusive-bucket-signal-contrast-frozen"
  )
  expect_identical(
    contrast$numeric_status,
    "threshold-only-vs-residual-only-signal-se-contrast-machine-readable"
  )
  expect_identical(
    contrast$consumer_summary_source,
    "story_local_wcb_partition_exclusive_bucket_representatives()"
  )
  expect_identical(contrast$threshold_only_replay_seeds, c(72L, 74L))
  expect_identical(contrast$residual_only_replay_seeds, c(48L, 53L, 79L))

  expect_identical(
    contrast$threshold_only_pair$closest_threshold_only$replay_seed,
    72L
  )
  expect_identical(
    contrast$threshold_only_pair$closest_threshold_only$treated_cluster_count,
    8L
  )
  expect_identical(
    contrast$threshold_only_pair$highest_t_threshold_only$replay_seed,
    74L
  )
  expect_identical(
    contrast$threshold_only_pair$highest_t_threshold_only$treated_cluster_count,
    9L
  )
  expect_identical(
    contrast$residual_only_pair$easiest_residual_only$replay_seed,
    48L
  )
  expect_identical(
    contrast$residual_only_pair$hardest_residual_only$replay_seed,
    79L
  )
  expect_identical(
    contrast$residual_only_pair$easiest_residual_only$treated_cluster_count,
    7L
  )
  expect_identical(
    contrast$residual_only_pair$hardest_residual_only$treated_cluster_count,
    7L
  )

  expect_equal(
    contrast$summary$bucket9_abs_att_gain_over_bucket8,
    2.175830683187882,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket9_original_se_reduction_vs_bucket8,
    0.602604466872937,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket9_t_stat_gain_over_bucket8,
    4.590259270017282,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket9_pvalue_gap_penalty_vs_bucket8,
    0.1875,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket8_abs_att_gain_over_easiest_residual,
    1.995846474090279,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket8_original_se_penalty_vs_easiest_residual,
    0.690099219619572,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket8_t_stat_gain_over_easiest_residual,
    1.110600598370500,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket9_abs_att_gain_over_hardest_residual,
    4.576686354486238,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket9_original_se_advantage_over_hardest_residual,
    0.401942677247878,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket9_t_stat_gain_over_hardest_residual,
    6.394904968579555,
    tolerance = 1e-15
  )
  expect_true(isTRUE(
    contrast$summary$all_threshold_only_abs_att_exceed_residual_only_family
  ))
  expect_true(isTRUE(
    contrast$summary$all_threshold_only_t_stats_exceed_residual_only_family
  ))
  expect_true(isTRUE(
    contrast$summary$bucket9_pairs_more_signal_with_lower_se_than_bucket8
  ))
  expect_true(isTRUE(
    contrast$summary$bucket8_retains_higher_t_than_easiest_residual_despite_se_penalty
  ))
  expect_true(isTRUE(
    contrast$summary$bucket9_dominates_hardest_residual_on_signal_and_precision
  ))
})

test_that("story-local TC-9.4.19 residual-only monotone decay keeps the 48/53/79 ladder machine-readable", {
  decay <- lwdid:::.story_local_wcb_partition_exclusive_residual_monotone_decay()

  expect_identical(decay$case_id, "TC-9.4.19")
  expect_identical(
    decay$exact_status,
    "story-local-partition-exclusive-residual-monotone-decay-frozen"
  )
  expect_identical(
    decay$numeric_status,
    "residual-only-48-53-79-monotone-decay-machine-readable"
  )
  expect_identical(
    decay$consumer_summary_source,
    "story_local_wcb_partition_exclusive_bucket_signal_contrast()"
  )
  expect_identical(decay$threshold_only_replay_seeds, c(72L, 74L))
  expect_identical(decay$residual_only_replay_seeds, c(48L, 53L, 79L))
  expect_identical(
    decay$threshold_only_anchor_pair$closest_threshold_only$replay_seed,
    72L
  )
  expect_identical(
    decay$threshold_only_anchor_pair$highest_t_threshold_only$replay_seed,
    74L
  )
  expect_identical(
    decay$residual_only_triplet$easiest_residual_only$replay_seed,
    48L
  )
  expect_identical(
    decay$residual_only_triplet$middle_residual_only$replay_seed,
    53L
  )
  expect_identical(
    decay$residual_only_triplet$hardest_residual_only$replay_seed,
    79L
  )
  expect_identical(
    unname(vapply(decay$residual_only_triplet, `[[`, integer(1L), "case_rank")),
    c(22L, 30L, 31L)
  )
  expect_identical(
    unname(vapply(decay$residual_only_triplet, `[[`, integer(1L), "replay_seed_offset")),
    c(6L, 11L, 37L)
  )
  expect_identical(
    unname(vapply(decay$residual_only_triplet, `[[`, integer(1L), "replay_attempt_id")),
    c(7L, 12L, 38L)
  )

  expect_equal(
    decay$summary$bucket48_abs_att_gain_over_bucket53,
    0.338242182004025,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket53_abs_att_gain_over_bucket79,
    0.066767015204051,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket53_original_se_penalty_vs_bucket48,
    0.276974687492090,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket53_abs_att_retention_vs_bucket48,
    0.299663567352167,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket53_original_se_multiplier_vs_bucket48,
    1.437870272295620,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket79_original_se_penalty_vs_bucket53,
    0.212462742502423,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket79_abs_att_retention_vs_bucket53,
    0.538675022988613,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket79_original_se_multiplier_vs_bucket53,
    1.233597644110348,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket48_t_stat_gain_over_bucket53,
    0.604404741423303,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket53_t_stat_retention_vs_bucket48,
    0.208407930204818,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket53_t_stat_gain_over_bucket79,
    0.089640358768470,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket79_t_stat_retention_vs_bucket53,
    0.436669951146914,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket53_pvalue_gap_penalty_vs_bucket48,
    0.4052734375,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket53_pvalue_gap_multiplier_vs_bucket48,
    1.924688057040998,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket79_pvalue_gap_penalty_vs_bucket53,
    0.0556640625,
    tolerance = 1e-12
  )
  expect_equal(
    decay$summary$bucket79_pvalue_gap_multiplier_vs_bucket53,
    1.065987497105812,
    tolerance = 1e-12
  )
  expect_true(isTRUE(decay$summary$residual_only_abs_att_strictly_descends))
  expect_true(isTRUE(decay$summary$residual_only_t_stats_strictly_descend))
  expect_true(isTRUE(decay$summary$residual_only_original_se_strictly_ascends))
  expect_true(isTRUE(decay$summary$residual_only_gap_above_threshold_strictly_ascends))
  expect_true(isTRUE(decay$summary$case_ranks_strictly_ascend))
  expect_true(isTRUE(decay$summary$replay_seed_offsets_strictly_ascend))
  expect_true(isTRUE(decay$summary$replay_attempt_ids_strictly_ascend))
  expect_true(isTRUE(
    decay$summary$threshold_only_anchor_abs_att_still_exceeds_residual_family
  ))
  expect_true(isTRUE(
    decay$summary$threshold_only_anchor_t_stats_still_exceed_residual_family
  ))
})

test_that("story-local TC-9.4.19 residual-only targeted replay rank stability keeps the 48/53/79 triplet machine-readable", {
  rank_stability <-
    lwdid:::.story_local_wcb_partition_exclusive_residual_targeted_replay_rank_stability()

  expect_identical(rank_stability$case_id, "TC-9.4.19")
  expect_identical(
    rank_stability$exact_status,
    "story-local-partition-exclusive-residual-targeted-replay-rank-stability-frozen"
  )
  expect_identical(
    rank_stability$numeric_status,
    "residual-only-48-53-79-targeted-replay-rank-stability-machine-readable"
  )
  expect_identical(rank_stability$residual_only_replay_seeds, c(48L, 53L, 79L))
  expect_identical(rank_stability$replay_seed_offsets, c(6L, 11L, 37L))
  expect_identical(rank_stability$replay_attempt_ids, c(7L, 12L, 38L))
  expect_identical(rank_stability$replay_case_ranks, c(22L, 30L, 31L))
  expect_identical(rank_stability$treated_cluster_pattern, c(7L, 7L, 7L))
  expect_identical(rank_stability$helper_rejection_pattern, c(FALSE, FALSE, FALSE))
  expect_identical(rank_stability$weakest_miss$seed, 48L)
  expect_identical(rank_stability$weakest_miss$attempt_id, 7L)
  expect_identical(rank_stability$weakest_miss$replication_id, 7L)
  expect_identical(rank_stability$weakest_miss$replay_seed_offset, 6L)
  expect_identical(rank_stability$weakest_miss$treated_cluster_count, 7L)
  expect_identical(rank_stability$weakest_miss$reject, FALSE)
  expect_equal(rank_stability$weakest_miss$pvalue, 0.48828125, tolerance = 1e-12)
  expect_equal(
    rank_stability$weakest_miss$abs_gap_to_threshold,
    0.43828125,
    tolerance = 1e-12
  )
  expect_equal(
    rank_stability$triplet_rows,
    list(
      list(
        case_rank = 22L,
        replay_seed = 48L,
        replay_seed_offset = 6L,
        replay_attempt_id = 7L,
        replication_id = 7L,
        treated_cluster_count = 7L,
        reject = FALSE,
        pvalue = 0.48828125,
        abs_gap_to_threshold = 0.43828125,
        att = 0.482970992563102,
        original_se = 0.632549650013912,
        t_stat_original = 0.763530566418746
      ),
      list(
        case_rank = 30L,
        replay_seed = 53L,
        replay_seed_offset = 11L,
        replay_attempt_id = 12L,
        replication_id = 12L,
        treated_cluster_count = 7L,
        reject = FALSE,
        pvalue = 0.8935546875,
        abs_gap_to_threshold = 0.8435546875,
        att = 0.144728810559076,
        original_se = 0.909524337506003,
        t_stat_original = 0.159125824995443
      ),
      list(
        case_rank = 31L,
        replay_seed = 79L,
        replay_seed_offset = 37L,
        replay_attempt_id = 38L,
        replication_id = 38L,
        treated_cluster_count = 7L,
        reject = FALSE,
        pvalue = 0.94921875,
        abs_gap_to_threshold = 0.89921875,
        att = 0.0779617953550249,
        original_se = 1.12198708000843,
        t_stat_original = 0.0694854662269725
      )
    ),
    tolerance = 1e-12
  )
  expect_identical(
    vapply(rank_stability$triplet_rows, `[[`, integer(1L), "replication_id"),
    c(7L, 12L, 38L)
  )
  expect_identical(
    vapply(rank_stability$triplet_rows, `[[`, integer(1L), "treated_cluster_count"),
    c(7L, 7L, 7L)
  )
  expect_identical(
    vapply(rank_stability$triplet_rows, `[[`, logical(1L), "reject"),
    c(FALSE, FALSE, FALSE)
  )
  expect_equal(
    rank_stability$summary$abs_att_values,
    c(0.482970992563102, 0.144728810559076, 0.0779617953550249),
    tolerance = 1e-12
  )
  expect_equal(
    rank_stability$summary$original_se_values,
    c(0.632549650013912, 0.909524337506003, 1.12198708000843),
    tolerance = 1e-12
  )
  expect_equal(
    rank_stability$summary$t_stat_values,
    c(0.763530566418746, 0.159125824995443, 0.0694854662269725),
    tolerance = 1e-12
  )
  expect_equal(
    rank_stability$summary$pvalue_values,
    c(0.48828125, 0.8935546875, 0.94921875),
    tolerance = 1e-12
  )
  expect_equal(
    rank_stability$summary$pvalue_gap_values,
    c(0.43828125, 0.8435546875, 0.89921875),
    tolerance = 1e-12
  )
  expect_identical(rank_stability$summary$targeted_replay_count, 3L)
  expect_identical(rank_stability$summary$replay_seed_labels, c(48L, 53L, 79L))
  expect_identical(rank_stability$summary$replication_ids, c(7L, 12L, 38L))
  expect_true(isTRUE(rank_stability$summary$replay_attempt_ids_strictly_ascend))
  expect_true(isTRUE(rank_stability$summary$replay_seed_offsets_strictly_ascend))
  expect_true(isTRUE(rank_stability$summary$case_ranks_strictly_ascend))
  expect_true(isTRUE(rank_stability$summary$abs_att_strictly_descends))
  expect_true(isTRUE(rank_stability$summary$original_se_strictly_ascends))
  expect_true(isTRUE(rank_stability$summary$t_stat_strictly_descends))
  expect_true(isTRUE(rank_stability$summary$pvalue_strictly_ascends))
  expect_true(isTRUE(rank_stability$summary$pvalue_gap_strictly_ascends))
  expect_true(isTRUE(rank_stability$summary$all_targeted_replays_remain_non_rejections))
  expect_true(isTRUE(rank_stability$summary$triplet_matches_residual_monotone_helper))
})

test_that("story-local TC-9.4.19 diagnostics reject invalid hardening thresholds", {
  threshold_error <- "`threshold` must be a single finite number in [0, 1)."

  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(threshold = NA_real_),
    threshold_error,
    fixed = TRUE
  )
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(threshold = 1),
    threshold_error,
    fixed = TRUE
  )
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(threshold = Inf),
    threshold_error,
    fixed = TRUE
  )
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(threshold = c(0.5, 0.6)),
    threshold_error,
    fixed = TRUE
  )
})

test_that("story-local TC-9.4.19 diagnostics reject unsupported weight distributions", {
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(weight_type = "gaussian"),
    paste(
      "`weight_type` must be one of",
      "'rademacher', 'mammen', or 'webb'."
    ),
    fixed = TRUE
  )
})

test_that("story-local TC-9.4.19 threshold validation makes the upper-exclusive bound explicit", {
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(threshold = 1),
    "`threshold` must be a single finite number in [0, 1).",
    fixed = TRUE
  )
})

test_that("story-local TC-9.4.19 diagnostics reject non-positive simulation counts", {
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(n_simulations = 0L),
    "`n_simulations` must be a single positive integer"
  )
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(n_simulations = -1L),
    "`n_simulations` must be a single positive integer"
  )
})

test_that("story-local TC-9.4.19 diagnostics reject invalid seed values", {
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(seed = -1L),
    "`seed` must be a single non-negative integer"
  )
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(seed = 1.5),
    "`seed` must be a single non-negative integer"
  )
})

test_that("story-local TC-9.4.19 diagnostics reject replay seeds before the canonical base seed", {
  expect_error(
    lwdid:::.story_local_wcb_power_diagnostics(replay_seeds = c(41L, 42L, 88L)),
    paste(
      "`replay_seeds` must map to canonical sequential attempts",
      "at or after the scenario base seed."
    )
  )
})

test_that("TC-9.4.23: bootstrap SE diagnostics use the population-SD formula on finite draws", {
  fixture <- make_wcb_fixture()

  result <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 123L,
      impose_null = TRUE,
      use_fwildclusterboot = FALSE
    )
  )
  att_valid <- result$att_bootstrap[is.finite(result$att_bootstrap)]
  se_diagnostics <- lwdid:::.compute_bootstrap_se_diagnostics(att_valid)

  expect_s3_class(result, "lwdid_wcb_result")
  expect_equal(length(att_valid), result$actual_n_bootstrap)
  expect_true(is.finite(se_diagnostics$se_bootstrap))
  expect_equal(result$se_bootstrap, se_diagnostics$se_bootstrap, tolerance = 1e-12)
  expect_equal(
    result$se_bootstrap,
    sqrt(mean((att_valid - mean(att_valid))^2)),
    tolerance = 1e-12
  )
})

test_that("TC-9.4.21: G=12 auto full enumeration is deterministic across seeds and requests", {
  fixture <- make_wcb_g12_fixture()

  run_123 <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 199L,
      weight_type = "rademacher",
      seed = 123L,
      use_fwildclusterboot = FALSE
    )
  )
  run_999 <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 999L,
      weight_type = "rademacher",
      seed = 999L,
      use_fwildclusterboot = FALSE
    )
  )
  run_forced <- local_suppress_lwdid_warnings(
    lwdid::wild_cluster_bootstrap(
      data = fixture,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 199L,
      weight_type = "rademacher",
      seed = 42L,
      full_enumeration = TRUE,
      use_fwildclusterboot = FALSE
    )
  )

  for (result in list(run_123, run_999, run_forced)) {
    expect_s3_class(result, "lwdid_wcb_result")
    expect_identical(result$weight_type, "rademacher")
    expect_identical(result$method, "native")
    expect_true(isTRUE(result$full_enumeration))
    expect_identical(result$actual_n_bootstrap, 4096L)
    expect_false(isTRUE(result$use_fwildclusterboot))
  }

  expect_equal(run_123$att, run_999$att, tolerance = 1e-12)
  expect_equal(run_123$att, run_forced$att, tolerance = 1e-12)
  expect_equal(run_123$original_se, run_999$original_se, tolerance = 1e-12)
  expect_equal(run_123$original_se, run_forced$original_se, tolerance = 1e-12)
  expect_equal(run_123$se_bootstrap, run_999$se_bootstrap, tolerance = 1e-12)
  expect_equal(run_123$se_bootstrap, run_forced$se_bootstrap, tolerance = 1e-12)
  expect_equal(run_123$pvalue, run_999$pvalue, tolerance = 1e-12)
  expect_equal(run_123$pvalue, run_forced$pvalue, tolerance = 1e-12)
  expect_equal(run_123$ci_lower, run_999$ci_lower, tolerance = 1e-12)
  expect_equal(run_123$ci_lower, run_forced$ci_lower, tolerance = 1e-12)
  expect_equal(run_123$ci_upper, run_999$ci_upper, tolerance = 1e-12)
  expect_equal(run_123$ci_upper, run_forced$ci_upper, tolerance = 1e-12)
  expect_equal(run_123$att, 1.9616666666666669, tolerance = 1e-12)
  expect_lt(abs(run_123$original_se - 0.01625151393834295), 1e-12)
  expect_equal(run_123$se_bootstrap, 0.5664920073967896, tolerance = 1e-12)
  expect_equal(run_123$pvalue, 0.000244140625, tolerance = 1e-12)
  expect_equal(run_123$ci_lower, 1.9300972764685163, tolerance = 1e-12)
  expect_equal(run_123$ci_upper, 1.9932360568648175, tolerance = 1e-12)
})

test_that("TC-9.4.26: requesting fwildclusterboot warns and falls back to the native path", {
  fixture <- make_wcb_fixture()

  result <- NULL
  expect_warning(
    result <- local_suppress_lwdid_warnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "Y",
        d = "D",
        cluster_var = "cluster",
        n_bootstrap = 99L,
        weight_type = "rademacher",
        seed = 1L,
        use_fwildclusterboot = TRUE
      )
    ),
    regexp = "fwildclusterboot",
    class = "lwdid_data"
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_equal(result$method, "native")
  expect_false(isTRUE(result$use_fwildclusterboot))
})
