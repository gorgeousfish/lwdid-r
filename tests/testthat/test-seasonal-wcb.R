make_seasonal_wcb_shared_dgp <- function() {
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

make_seasonal_wcb_decision_scenario <- function(
    scenario_id,
    true_tau,
    n_simulations,
    metric_name,
    decision_direction,
    N = 200L,
    G = 10L,
    seed = 42L
) {
  if (N %% G != 0L) {
    stop("N must be divisible by G for the seasonal WCB fixture scaffold.", call. = FALSE)
  }

  list(
    id = scenario_id,
    seed = as.integer(seed),
    n_simulations = as.integer(n_simulations),
    true_tau = as.numeric(true_tau),
    G = as.integer(G),
    obs_per_cluster = as.integer(N / G),
    estimator = "wild_cluster_bootstrap",
    n_bootstrap = 199L,
    requested_n_bootstrap = 199L,
    weight_type = "rademacher",
    alpha = 0.05,
    impose_null = TRUE,
    threshold = 0.70,
    metric_name = metric_name,
    decision_direction = decision_direction,
    pvalue_threshold = 0.05
  )
}

simulate_clustered_effect_fixture <- function(
    N = 200L,
    G = 10L,
    true_tau,
    seed = 42L,
    scenario_id,
    metric_name,
    decision_direction
) {
  scenario <- make_seasonal_wcb_decision_scenario(
    scenario_id = scenario_id,
    true_tau = true_tau,
    n_simulations = 1L,
    metric_name = metric_name,
    decision_direction = decision_direction,
    N = N,
    G = G,
    seed = seed
  )

  set.seed(as.integer(seed))
  lwdid:::.simulate_clustering_monte_carlo_panel(
    scenario = scenario,
    shared_dgp = make_seasonal_wcb_shared_dgp()
  )
}

simulate_zero_effect_clustered <- function(N = 200L, G = 10L, seed = 42L) {
  simulate_clustered_effect_fixture(
    N = N,
    G = G,
    true_tau = 0,
    seed = seed,
    scenario_id = "story_e9_06_zero_effect_fixture",
    metric_name = "non_rejection_rate",
    decision_direction = "gt"
  )
}

simulate_strong_effect_clustered <- function(
    N = 200L,
    G = 10L,
    tau = 2.0,
    seed = 42L
) {
  simulate_clustered_effect_fixture(
    N = N,
    G = G,
    true_tau = tau,
    seed = seed,
    scenario_id = "story_e9_06_strong_effect_fixture",
    metric_name = "power_under_alternative",
    decision_direction = "lt"
  )
}

run_seasonal_wcb_public_smoke <- function(seed = 42L) {
  strong_effect <- simulate_strong_effect_clustered(seed = seed)

  lwdid::wild_cluster_bootstrap(
    data = strong_effect,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    n_bootstrap = 199L,
    weight_type = "rademacher",
    seed = as.integer(seed),
    use_fwildclusterboot = FALSE
  )
}

make_e906_quarterly_panel <- function() {
  data.frame(
    id = rep(1L, 6L),
    tindex = 1:6,
    y = c(10, 11, 12, 13, 14, 15),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L)
  )
}

make_e906_monthly_panel <- function() {
  data.frame(
    id = rep(1L, 14L),
    tindex = 1:14,
    y = 100 + seq_len(14L),
    month = c(1:12, 1L, 2L),
    post = c(rep(0L, 13L), 1L)
  )
}

make_e906_weekly_panel <- function() {
  data.frame(
    id = rep(1L, 54L),
    tindex = 1:54,
    y = 200 + seq_len(54L),
    week = c(1:52, 1L, 2L),
    post = c(rep(0L, 53L), 1L)
  )
}

make_e906_weekly_detrend_panel <- function() {
  week <- c(1:52, 1L, 2L, 3L)
  tindex <- seq_along(week)

  data.frame(
    id = rep(1L, length(week)),
    tindex = tindex,
    week = week,
    post = c(rep(0L, 54L), 1L),
    y = 50 + 0.2 * tindex + sin(2 * pi * week / 52)
  )
}

make_e906_frequency_detection_cases <- function() {
  list(
    quarterly = list(
      data = data.frame(
        id = rep(1:2, each = 8L),
        year = rep(rep(2001:2002, each = 4L), times = 2L),
        quarter = rep(rep(1:4, times = 2L), times = 2L),
        y = seq_len(16L)
      ),
      expected_frequency = "quarterly",
      expected_Q = 4L
    ),
    monthly = list(
      data = data.frame(
        id = rep(1:2, each = 24L),
        year = rep(rep(2001:2002, each = 12L), times = 2L),
        month = rep(rep(1:12, times = 2L), times = 2L),
        y = seq_len(48L)
      ),
      expected_frequency = "monthly",
      expected_Q = 12L
    ),
    weekly = list(
      data = data.frame(
        id = rep(1:2, each = 104L),
        year = rep(rep(2001:2002, each = 52L), times = 2L),
        week = rep(rep(1:52, times = 2L), times = 2L),
        y = seq_len(208L)
      ),
      expected_frequency = "weekly",
      expected_Q = 52L
    )
  )
}

make_e906_validation_case_bundle <- function() {
  invalid_range_panel <- make_e906_quarterly_panel()
  invalid_range_panel$quarter[1] <- 5L

  min_global_pre_panel <- data.frame(
    id = rep(1L, 3L),
    tindex = 1:3,
    quarter = c(1L, 2L, 3L),
    post = c(0L, 1L, 1L),
    y = c(10.0, 11.0, 12.0)
  )

  list(
    a1_invalid_range = list(
      call = list(
        data = invalid_range_panel,
        y = "y",
        season_var = "quarter",
        Q = 4L,
        ivar = "id",
        tindex = "tindex",
        post = "post",
        method = "demeanq",
        exclude_pre_periods = 0L,
        min_global_pre_periods = 1L
      ),
      regexp = "contains invalid values|Expected integer values"
    ),
    a2_coverage_warning = list(
      call = list(
        data = make_e906_missing_q4_panel(),
        y = "y",
        season_var = "quarter",
        Q = 4L,
        ivar = "id",
        tindex = "tindex",
        post = "post",
        method = "demeanq",
        exclude_pre_periods = 0L,
        min_global_pre_periods = 1L
      ),
      regexp = "not observed in the pre-treatment window|seasonal gate to pass"
    ),
    a3_insufficient_pre = list(
      call = list(
        data = make_e906_sparse_panel(),
        y = "y",
        season_var = "quarter",
        Q = 4L,
        ivar = "id",
        tindex = "tindex",
        post = "post",
        method = "demeanq",
        exclude_pre_periods = 0L,
        min_global_pre_periods = 1L
      ),
      regexp = "requires at least 3 observation\\(s\\)|df >= 1"
    ),
    a6_exclude_pre_base = list(
      call = list(
        data = make_e906_quarterly_panel(),
        y = "y",
        season_var = "quarter",
        Q = 4L,
        ivar = "id",
        tindex = "tindex",
        post = "post",
        method = "demeanq",
        exclude_pre_periods = 0L,
        min_global_pre_periods = 1L
      )
    ),
    a6_exclude_pre_shifted = list(
      call = list(
        data = make_e906_quarterly_panel(),
        y = "y",
        season_var = "quarter",
        Q = 4L,
        ivar = "id",
        tindex = "tindex",
        post = "post",
        method = "demeanq",
        exclude_pre_periods = 2L,
        min_global_pre_periods = 1L
      ),
      regexp = "requires at least 4 observation\\(s\\)|df >= 1"
    ),
    a7_min_global_pre = list(
      call = list(
        data = min_global_pre_panel,
        y = "y",
        season_var = "quarter",
        Q = 4L,
        ivar = "id",
        tindex = "tindex",
        post = "post",
        method = "detrendq",
        exclude_pre_periods = 0L,
        min_global_pre_periods = 2L
      ),
      regexp = "requires at least 2 pre-treatment period\\(s\\)|Found: 1"
    ),
    a8_quarter_alias = list(
      call = list(
        data = make_e906_quarterly_panel(),
        y = "y",
        season_var = NULL,
        Q = 4L,
        ivar = "id",
        tindex = "tindex",
        post = "post",
        method = "demeanq",
        exclude_pre_periods = 0L,
        min_global_pre_periods = 1L,
        quarter = "quarter"
      )
    )
  )
}

make_e906_known_seasonal_panel <- function() {
  quarter <- c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L)
  epsilon <- c(rep(0.0, 8L), 0.10, -0.20, 0.30, -0.40)
  expected_seasonal_fit <- unname(c(`1` = 0, `2` = 2, `3` = 3, `4` = 4)[as.character(quarter)])

  data.frame(
    id = rep(1L, length(quarter)),
    tindex = seq_along(quarter),
    quarter = quarter,
    post = c(rep(0L, 8L), rep(1L, 4L)),
    epsilon = epsilon,
    expected_seasonal_fit = expected_seasonal_fit,
    y = 10 + expected_seasonal_fit + epsilon
  )
}

make_e906_known_trend_seasonal_panel <- function() {
  quarter <- c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L)
  tindex <- seq_along(quarter)
  epsilon <- c(0.00, 0.10, -0.10, 0.20, -0.05, 0.15, -0.15, 0.05, 0.00, 0.10)
  expected_trend_component <- 0.3 * tindex
  expected_seasonal_fit <- unname(c(`1` = 0, `2` = 2, `3` = 3, `4` = 4)[as.character(quarter)])

  data.frame(
    id = rep(1L, length(quarter)),
    tindex = tindex,
    quarter = quarter,
    post = c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    epsilon = epsilon,
    expected_trend_component = expected_trend_component,
    expected_seasonal_fit = expected_seasonal_fit,
    y = 5 + expected_trend_component + expected_seasonal_fit + epsilon
  )
}

make_e906_missing_q4_panel <- function() {
  data.frame(
    id = rep(1L, 7L),
    tindex = seq_len(7L),
    quarter = c(1L, 2L, 1L, 2L, 1L, 3L, 4L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L, 1L),
    y = c(10, 12, 10, 12, 10, 16, 17)
  )
}

make_e906_single_cluster_panel <- function() {
  data.frame(
    Y = c(5.0, 5.5, 6.0, 6.5, 7.0, 7.5),
    D = c(0L, 0L, 0L, 1L, 1L, 1L),
    cluster = rep(1L, 6L),
    post = c(0L, 0L, 0L, 1L, 1L, 1L)
  )
}

make_e906_constant_outcome_panel <- function() {
  data.frame(
    Y = rep(7.0, 12L),
    D = rep(c(0L, 1L), each = 6L),
    cluster = rep(1:3, each = 4L),
    post = rep(c(0L, 0L, 1L, 1L), times = 3L)
  )
}

make_e906_constant_treatment_cluster_panel <- function() {
  cluster <- rep(1:4, each = 2L)
  control_a <- rep(c(0, 1), times = 4L)
  control_b <- rep(c(-1, 1), times = 4L)
  cluster_effect <- c(-0.4, 0.2, -0.1, 0.3)[cluster]
  within_effect <- rep(c(-0.2, 0.15), times = 4L)

  data.frame(
    Y = 2.5 + 0.5 * control_a - 0.35 * control_b + cluster_effect + within_effect,
    D = rep(1L, length(cluster)),
    cluster = cluster,
    post = rep(c(0L, 1L), times = 4L),
    control_a = control_a,
    control_b = control_b
  )
}

make_e906_sparse_panel <- function() {
  data.frame(
    id = rep(1L, 3L),
    tindex = 1:3,
    quarter = c(1L, 2L, 3L),
    post = c(0L, 0L, 1L),
    y = c(10.0, 10.5, 11.0)
  )
}

make_e906_constant_time_panel <- function() {
  data.frame(
    id = rep(1:2, each = 7L),
    tindex = c(
      rep(1L, 5L), 6L, 7L,
      rep(2L, 5L), 6L, 7L
    ),
    quarter = rep(c(1L, 2L, 3L, 1L, 2L, 3L, 4L), times = 2L),
    post = rep(c(0L, 0L, 0L, 0L, 0L, 1L, 1L), times = 2L),
    y = c(
      6.0, 6.5, 7.0, 6.2, 6.7, 8.0, 8.5,
      9.0, 9.5, 10.0, 9.2, 9.7, 11.0, 11.5
    )
  )
}

make_e906_all_sparse_panel <- function() {
  data.frame(
    id = rep(1:2, each = 3L),
    tindex = rep(1:3, times = 2L),
    quarter = c(1L, 2L, 3L, 1L, 2L, 4L),
    post = rep(c(0L, 0L, 1L), times = 2L),
    y = c(9.0, 9.5, 10.0, 11.0, 11.5, 12.0)
  )
}

make_e906_biannual_panel <- function() {
  data.frame(
    id = rep(1L, 6L),
    tindex = 1:6,
    halfyear = c(1L, 2L, 1L, 2L, 1L, 2L),
    post = c(0L, 0L, 0L, 0L, 1L, 1L),
    y = c(20.0, 21.0, 20.5, 21.5, 22.0, 23.0)
  )
}

make_e906_biannual_identity_panel <- function() {
  data.frame(
    id = rep(1:2, each = 10L),
    tindex = rep(1:10, times = 2L),
    halfyear = rep(rep(c(1L, 2L), times = 5L), times = 2L),
    post = rep(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L), times = 2L),
    y = c(
      20.0, 21.0, 20.5, 21.5, 20.25, 21.25, 22.0, 23.0, 22.5, 23.5,
      30.0, 31.0, 30.5, 31.5, 30.25, 31.25, 32.0, 33.0, 32.5, 33.5
    ),
    expected_seasonal_fit = c(
      rep(c(20.25, 21.25), times = 5L),
      rep(c(30.25, 31.25), times = 5L)
    ),
    expected_ydot = c(
      -0.25, -0.25, 0.25, 0.25, 0.00, 0.00, 1.75, 1.75, 2.25, 2.25,
      -0.25, -0.25, 0.25, 0.25, 0.00, 0.00, 1.75, 1.75, 2.25, 2.25
    )
  )
}

make_e906_small_ti_panel <- function() {
  cluster <- rep(1:3, each = 4L)
  post <- rep(c(0L, 0L, 1L, 1L), times = 3L)
  D <- as.integer(cluster >= 2L & post == 1L)

  data.frame(
    Y = 4.0 + 0.15 * seq_len(12L) + 0.4 * cluster + 1.25 * D,
    D = D,
    cluster = cluster,
    post = post
  )
}

make_e906_grid_resolution_panel <- function() {
  cluster_sizes <- c(13L, 13L, 13L, 13L, 12L, 12L, 12L, 12L)
  cluster <- rep(1:8, times = cluster_sizes)
  observation <- seq_along(cluster)
  D <- as.integer(rep(c(0L, 1L, 0L, 1L, 1L), length.out = length(cluster)))
  cluster_effect <- c(-0.45, 0.12, -0.30, 0.28, -0.18, 0.34, -0.09, 0.21)[cluster]
  epsilon <- 0.6 * sin(observation / 3) + 0.35 * cos(observation / 5)

  data.frame(
    Y = 1.0 + 2.0 * D + cluster_effect + epsilon,
    D = D,
    cluster = cluster
  )
}

make_e906_core_wcb_panel <- function() {
  cluster <- rep(1:8, each = 4L)
  obs_within_cluster <- rep(1:4, times = 8L)
  D <- as.integer(cluster %in% c(2L, 4L, 6L, 8L))
  cluster_effect <- c(-0.45, 0.20, -0.35, 0.40, -0.15, 0.30, -0.25, 0.35)[cluster]
  within_effect <- c(-0.18, 0.04, 0.12, -0.02)[obs_within_cluster]

  data.frame(
    Y = 1.0 + 0.8 * D + cluster_effect + within_effect + 0.01 * cluster * obs_within_cluster,
    D = D,
    cluster = cluster
  )
}

make_e906_public_full_enumeration_panel <- function() {
  make_e906_core_wcb_panel()
}

make_e906_public_g12_full_enumeration_panel <- function() {
  cluster <- rep(1:12, each = 3L)
  D <- rep(c(0L, 0L, 1L), times = 12L)
  cluster_level <- rep(seq(-1.1, 1.1, length.out = 12L), each = 3L)
  within_cluster <- rep(c(-0.25, 0.15, 0.35), times = 12L)
  deterministic_noise <- c(
    rep(c(-0.10, 0.05, 0.20), 4L),
    rep(c(0.08, -0.03, 0.12), 4L),
    rep(c(-0.06, 0.09, 0.18), 4L)
  )

  data.frame(
    Y = 2.0 + 0.6 * cluster_level + 1.4 * D + within_cluster + deterministic_noise,
    D = D,
    cluster = cluster
  )
}

make_e906_panel_data <- function() {
  id <- rep(1:4, each = 10L)
  tindex <- rep(1:10, times = 4L)
  quarter <- rep(c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L), times = 4L)
  post <- as.integer(tindex >= 9L)
  D <- as.integer(id >= 3L & post == 1L)
  unit_effect <- c(0.0, 0.2, 0.5, 0.7)[id]
  seasonal_effect <- c(0.0, 0.4, 0.8, 1.2)[quarter]

  data.frame(
    id = id,
    tindex = tindex,
    quarter = quarter,
    post = post,
    D = D,
    y = 5.0 + unit_effect + 0.3 * tindex + seasonal_effect + 1.5 * D
  )
}

make_e906_collinear_panel <- function() {
  panel_data <- make_e906_panel_data()
  panel_data$control_a <- seq_len(nrow(panel_data))
  panel_data$control_b <- 2 * panel_data$control_a
  panel_data
}

make_e906_restricted_model_panel <- function() {
  data.frame(
    cluster = rep(1:8, each = 4L),
    D = c(
      0L, 0L, 0L, 0L,
      1L, 1L, 1L, 1L,
      0L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L,
      1L, 1L, 1L, 1L,
      0L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L
    ),
    control_a = c(
      0.32950777, -0.82046838, 0.48742905, 0.73832471,
      0.57578135, -0.30538839, 1.51178117, 0.38984324,
      -0.62124058, -2.21469989, 1.12493092, -0.04493361,
      -0.01619026, 0.94383621, 0.82122120, 0.59390132,
      0.91897737, 0.78213630, 0.07456498, -1.98935170,
      0.61982575, -0.05612874, -0.15579551, -1.47075238,
      -0.47815006, 0.41794156, 1.35867955, -0.10278773,
      0.38767161, -0.05380504, -1.37705956, -0.41499456
    ),
    control_b = c(
      -0.20778069, -0.57859065, 0.09893121, 1.00589408,
      0.90141335, -0.31163743, 1.00817575, 1.03818897,
      -0.31705320, -0.23938061, 0.13746801, -0.36277578,
      0.50931722, -0.31146463, 0.66860880, -0.57844576,
      -0.23332976, -0.47875035, -0.49837167, -1.47893675,
      0.47052424, 0.73569980, 0.51109070, 0.15339194,
      -0.23289611, -0.05444937, 1.02934435, 0.17903026,
      0.42574934, -0.30974697, -0.87259758, 0.86086975
    ),
    y = c(
      5.139622, 3.909829, 4.668841, 4.225553,
      4.288904, 4.511627, 4.625301, 4.642567,
      4.621786, 3.972368, 5.697409, 4.832262,
      5.074930, 5.771133, 5.115382, 5.966055,
      5.628793, 5.766814, 5.299634, 4.179209,
      3.429968, 2.932585, 3.303439, 1.886405,
      4.272061, 4.744853, 5.044010, 4.104319,
      4.069967, 4.102952, 3.255694, 3.500003
    )
  )
}

make_e906_test_inversion_panel <- function() {
  make_e906_restricted_model_panel()[, c("y", "D", "cluster")]
}

make_e906_seasonal_data <- function() {
  make_e906_panel_data()
}

make_e906_clustered_data <- function() {
  clustered_data <- make_e906_panel_data()
  clustered_data$cluster <- rep(
    c(1L, 1L, 2L, 2L),
    each = max(tabulate(clustered_data$id))
  )
  clustered_data
}

make_e906_seasonal_clustered_data <- function() {
  make_e906_clustered_data()
}

make_e906_public_lwdid_panel <- function(with_anticipation = FALSE) {
  id <- rep(1:8, each = 12L)
  time <- rep(1:12, times = 8L)
  quarter <- rep(c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L), times = 8L)
  post <- as.integer(time >= 9L)
  treat <- as.integer(id >= 5L)
  cluster <- rep(rep(1:4, each = 2L), each = 12L)
  unit_effect <- c(-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0)[id]
  seasonal_effect <- c(0.0, 0.3, 0.6, 0.9)[quarter]
  id_noise <- c(0.05, -0.03, 0.02, -0.04, 0.01, -0.02, 0.03, -0.01)[id]
  time_noise <- c(0.00, 0.02, -0.01, 0.03, -0.02, 0.01, 0.04, -0.03, 0.05, -0.04, 0.02, -0.01)[time]
  anticipation_effect <- if (isTRUE(with_anticipation)) {
    as.numeric(treat == 1L & time %in% c(7L, 8L)) * 0.8
  } else {
    0.0
  }

  data.frame(
    id = id,
    time = time,
    quarter = quarter,
    post = post,
    treat = treat,
    cluster = cluster,
    y = 5.0 +
      unit_effect +
      0.25 * time +
      seasonal_effect +
      1.25 * (treat * post) +
      anticipation_effect +
      id_noise * time_noise * 10
  )
}

make_e906_public_lwdid_anticipation_panel <- function() {
  make_e906_public_lwdid_panel(with_anticipation = TRUE)
}

make_e906_post_consume_closure_wave_public_contracts <- function() {
  list(
    quarter_alias_demeanq_wcb = list(
      att = 1.25,
      pvalue = 0.125,
      ci_lower = -0.01940512758846036,
      ci_upper = 2.519405127588462,
      actual_n_bootstrap = 16L,
      requested_n_bootstrap = 19L
    ),
    detrendq_wcb = list(
      att = 1.250375,
      pvalue = 0.125,
      ci_lower = 0,
      ci_upper = 2.500750000000001,
      actual_n_bootstrap = 16L,
      requested_n_bootstrap = 19L
    ),
    detrendq_bootstrap_anticipation = list(
      att = 0.4503749999999988,
      pvalue = 0.0625,
      ci_lower = 0,
      ci_upper = 0.9007499999999976,
      actual_n_bootstrap = 16L,
      requested_n_bootstrap = 19L
    ),
    detrendq_excluded_auto_wcb = list(
      att = 1.250624999999999,
      pvalue = 0.125,
      ci_lower = -2.220446049250313e-16,
      ci_upper = 2.501249999999998,
      actual_n_bootstrap = 16L,
      requested_n_bootstrap = 19L
    )
  )
}

make_e906_layer3_realdata_contracts <- function() {
  list(
    quarterly = jsonlite::fromJSON(
      "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e9_01_layer3_smoking_quarterly_realdata_regression.json",
      simplifyVector = FALSE
    ),
    small_g = jsonlite::fromJSON(
      "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e9_04_layer3_smoking_small_g_public_replay.json",
      simplifyVector = FALSE
    )
  )
}

make_e906_smoking_quarterly_realdata_panel <- function() {
  utils::read.csv(
    "/Users/cxy/Desktop/lwdid_r/lwdid-py_v0.2.3/tests/data/smoking_quarterly.csv",
    stringsAsFactors = FALSE
  )
}

make_e906_smoking_small_g_public_wcb_fixture <- function() {
  utils::data("smoking", package = "lwdid", envir = environment())

  state_subset <- sort(unique(smoking$state))[1:10]
  smoking_small_g <- subset(smoking, state %in% state_subset)
  transformed <- suppressWarnings(
    lwdid::lwdid(
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
    y_dot = transformed$.ri_data$y_trans,
    treated = transformed$.ri_data$d,
    state = state_subset
  )
}

make_e906_test_inversion_shared_fixture <- function() {
  utils::read.csv(
    "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity/e9_04_layer2_test_inversion_shared_controls_fixture.csv",
    stringsAsFactors = FALSE
  )
}

# Freeze the sampled six-scenario Python optimization reference for D10b inside
# the owner test bundle so the active file can exercise the comparator directly.
make_e906_python_multiscenario_wcb_fixture <- function() {
  data.frame(
    y_transformed = c(
      0.309128269584344, 1.933545479453588, 2.274206172948465, 2.818235099176404, 2.068261082233229, -0.174034177862348, 1.182990032471633, 0.668189106543637, 2.824685445515255, 4.155116365181820,
      3.690473976881142, 2.489756647319959, 2.997852860627315, 2.456415928206513, 2.928499722996716, 1.394855331148437, 1.328931898330927, 0.067970467643849, 1.242848321505090, 3.380586870016121,
      4.997481984413673, 4.356852558161422, -0.373235403679068, 0.192624585330576, -0.983434054957251, -0.110699067083734, 0.795612950470093, 1.685861567767986, -0.249075890863856, 1.825861506811955,
      -1.667281082858903, 2.317342026334845, 1.417593535020189, -0.049431480486322, 1.603147168712164, 2.385745392627195, -1.277720444476879, -0.987037298511238, 4.609255058908050, -0.807414634394279,
      0.278440019448589, 3.212139126224542, 6.280292767143635, 2.203658969077887, 3.006976641115198, 3.716781148843738, 3.103923242308821, 0.388368519474991, 1.129947751961119, 2.152563173964973,
      3.809991672369847, 0.999069526673129, 1.241402633387187, 2.000329789533247, 3.137051853928563, 1.108838107894274, 1.607395487529485, 0.819898379084213, 1.846764120448368, 2.387941415723057,
      1.625711685864073, 1.822038766594811, 3.523437705773474, 1.470366751287021, 0.748724324652760, 1.786776619868915, 1.372995759782118, 5.563450667126202, 3.043637260670850, 3.856128910752108,
      2.707225148500839, 0.357870886212189, 4.661425311284472, 1.733031583342746, 1.950425761587218, -0.274303628778586, 4.397720839753788, 1.924702696409970, 0.359851673781129, 0.998757550745658,
      0.288199610029601, 2.076949204029325, -0.267688780130473, -0.079005862392879, -0.227494972533804, 2.177044917775550, 1.495974378666239, -2.542176035258216, -0.212367858836602, 2.200994827045529,
      2.320476395360181, 1.652567824351792, -1.732239376805501, 2.352477889387191, 3.751832619837267, 0.490621000756944, 0.888310419371270, -0.304042817449812, -0.036841940155766, 0.971583894969133,
      0.014454074212964, 1.648954841267799, 0.931565488556548, 2.081221851978785, 1.526322135546535, 2.787996925498417, 4.092844772736633, 2.787903405307140, 1.176898423365900, 2.052339293034132,
      0.415456062428046, 1.651486423268564, 1.338494181163694, -0.829393376546025, 1.131260588986502, 0.158259189063681, 1.520157340036160, 1.099779235300101, 2.169556404215642, 1.079684937707486,
      -0.508120054069598, 2.996967028685442, -0.108086991011173, 0.545640359217368, 2.575148634860364, 1.701454012481339, 2.364989773004896, 0.050202835143244, 3.096702294459904, -0.152639947797485,
      0.392779924187062, 1.375929673134886, 4.219232237690338, 1.672681673808083, -1.227369307014457, 1.368286930613902, 1.135953774816221, 1.477810486956336, 2.680129883234850, -0.999318075221100,
      2.250775278637104, 3.432785026621564, -0.739715695625733, 2.317827751881718, 1.210235316528867, 3.688462689671613, 0.270935042306414, 3.655642393794420, 3.910350259817247, 3.074756723965377,
      1.076279596526905, 3.112221943109935, -0.021466560911173, 0.694350023587856, 0.689356171375850, 1.225272568181704, 3.841702950738504, 1.783342621849800, 0.716869673508631, 4.323334251487182,
      0.248011192435291, 0.169349403957535, 1.193905613079049, 3.835984611733176, 1.003025476808515, 2.320560011684845, 1.848393403814383, 3.833419488137917, 3.545715922571949, 5.319252521593940,
      2.786437732381811, 0.988077129181342, 2.595853782764997, 0.533899990291030, 3.099479865590889, 4.301150106108382, -0.610660119811070, -0.260047792858284, 1.390321484129002, 0.197463131222432,
      3.013114642230204, 1.097665879025981, 0.849752215541764, 1.160712152320843, -1.736625400419434, 0.627831621603592, 0.101702002574874, 3.032863807929822, 1.275431662434480, 0.711719468271892,
      3.333611121589064, 2.802689410328278, 1.103542741814677, -0.131195226656672, 0.393866622381751, 0.770254132870167, 1.394006511630629, 2.179736863988638, 1.094645291945322, 0.307639244297831
    ),
    d = c(
      0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
      1, 0, 1, 1, 1, 0, 0, 0, 0, 1,
      1, 1, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 1, 1, 0, 1, 1, 0, 0, 1, 0,
      0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
      1, 0, 0, 1, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 1, 0, 1,
      1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
      1, 1, 0, 1, 1, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 1, 1, 0, 1,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 1, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 1, 0, 0, 0, 1, 0, 1, 1, 1,
      0, 1, 0, 0, 0, 0, 1, 1, 0, 1,
      0, 0, 0, 1, 0, 0, 0, 1, 1, 1,
      1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
      1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
      1, 1, 0, 0, 0, 0, 0, 1, 0, 0
    ),
    cluster = rep(1:10, each = 20L)
  )
}

make_e906_python_multiscenario_wcb_scenarios <- function() {
  list(
    impose_null_true_rademacher = list(
      weight_type = "rademacher",
      impose_null = TRUE,
      att = 2.265783711291041,
      t_stat_original = 15.60019152571567,
      se_bootstrap = 0.7347627885478468,
      pvalue = 0.002002002002002002,
      ci_lower = 1.9219300686787044,
      ci_upper = 2.6096373539033775,
      n_bootstrap = 999L
    ),
    impose_null_true_mammen = list(
      weight_type = "mammen",
      impose_null = TRUE,
      att = 2.265783711291041,
      t_stat_original = 15.60019152571567,
      se_bootstrap = 0.6999231251956204,
      pvalue = 0.03003003003003003,
      ci_lower = 1.8326007164515576,
      ci_upper = 2.6989667061305247,
      n_bootstrap = 999L
    ),
    impose_null_true_webb = list(
      weight_type = "webb",
      impose_null = TRUE,
      att = 2.265783711291041,
      t_stat_original = 15.60019152571567,
      se_bootstrap = 0.7191167672205523,
      pvalue = 0,
      ci_lower = 1.9276286567743854,
      ci_upper = 2.6039387658076967,
      n_bootstrap = 999L
    ),
    impose_null_false_rademacher = list(
      weight_type = "rademacher",
      impose_null = FALSE,
      att = 2.265783711291041,
      t_stat_original = 15.60019152571567,
      se_bootstrap = 0.1396965308530579,
      pvalue = 0.6266266266266266,
      ci_lower = 2.0002998716806064,
      ci_upper = 2.5295300857655074,
      n_bootstrap = 999L
    ),
    impose_null_false_mammen = list(
      weight_type = "mammen",
      impose_null = FALSE,
      att = 2.265783711291041,
      t_stat_original = 15.60019152571567,
      se_bootstrap = 0.1315040216251721,
      pvalue = 0.6876876876876877,
      ci_lower = 2.0331379765235944,
      ci_upper = 2.535806835760164,
      n_bootstrap = 999L
    ),
    impose_null_false_webb = list(
      weight_type = "webb",
      impose_null = FALSE,
      att = 2.265783711291041,
      t_stat_original = 15.60019152571567,
      se_bootstrap = 0.1403637085345523,
      pvalue = 0.6106106106106106,
      ci_lower = 1.9916972151785377,
      ci_upper = 2.5335412267441964,
      n_bootstrap = 999L
    )
  )
}

# Freeze the Python-aligned oracle data inside the story-local suite so the
# active owner file can exercise B4/C4/D10 without depending on older parity files.
make_e906_python_reference_bundle <- function() {
  wcb_cluster_base <- c(1.2, 1.4, 0.9, 1.1, 1.3, 3.4, 3.7, 3.6, 3.9, 4.1)

  list(
    demeanq = list(
      source_anchor = "tests/test_transformations.py::TestQuarterlyTransforms.test_demeanq_mve",
      fixture = data.frame(
        id = rep(1L, 9L),
        year = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
        quarter = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L),
        y = c(10.0, 12.0, 11.0, 13.0, 10.0, 15.0, 17.0, 16.0, 18.0),
        d_ = rep(1L, 9L),
        post_ = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
        tindex = 1:9
      ),
      python_yhat = c(10.0, 12.0, 11.0, 13.0, 10.0, 12.0, 11.0, 13.0, 10.0),
      python_ydot = c(0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 6.0, 3.0, 8.0)
    ),
    detrendq = list(
      source_anchor = "tests/test_transformations.py::TestQuarterlyTransforms.test_detrendq_mve",
      fixture = data.frame(
        id = rep(1L, 10L),
        year = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
        quarter = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L),
        time = 1:10,
        post = c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
        y = c(5.5, 7.0, 7.0, 8.5, 7.5, 9.0, 9.0, 10.5, 10.0, 11.5)
      ),
      python_yhat = c(5.5, 7.0, 7.0, 8.5, 7.5, 9.0, 9.0, 10.5, 9.5, 11.0),
      python_ydot = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5)
    ),
    wild_bootstrap = list(
      source_anchor = "tests/inference/test_wild_bootstrap_optimization.py::_run_bootstrap",
      fixture = data.frame(
        y_transformed = unlist(
          lapply(wcb_cluster_base, function(mu) c(mu - 0.15, mu + 0.15)),
          use.names = FALSE
        ),
        d = c(rep(0, 10L), rep(1, 10L)),
        cluster = rep(1:10, each = 2L)
      ),
      python_att = 2.5600000000000005,
      python_t_stat_original = 17.818220876824675,
      python_pvalue = 0.001953125,
      python_n_bootstrap = 1024L,
      python_ci_lower = 2.2093772944268597,
      python_ci_upper = 2.9106227055731413,
      n_clusters = 10L
    ),
    wild_bootstrap_multiscenario = list(
      source_anchor = paste(
        "tests/inference/test_wild_bootstrap_optimization.py::SCENARIOS",
        "+ reference_data/wild_bootstrap_reference.json"
      ),
      fixture = make_e906_python_multiscenario_wcb_fixture(),
      n_bootstrap = 999L,
      alpha = 0.05,
      seed = 42L,
      n_clusters = 10L,
      shared_att = 2.265783711291041,
      shared_t_stat_original = 15.60019152571567,
      scenarios = make_e906_python_multiscenario_wcb_scenarios()
    )
  )
}

make_e906_python_transform_package_contracts <- function() {
  oracle_bundle <- make_e906_python_reference_bundle()

  list(
    demeanq = list(
      source_anchor = paste(
        oracle_bundle$demeanq$source_anchor,
        "+ transform_demeanq / transform_common package surface"
      ),
      seasonal_fit = oracle_bundle$demeanq$python_yhat,
      y_trans = oracle_bundle$demeanq$python_ydot,
      ydot_postavg = rep(5, length(oracle_bundle$demeanq$python_ydot)),
      firstpost = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
      n_pre = rep(5L, length(oracle_bundle$demeanq$python_ydot))
    ),
    detrendq = list(
      source_anchor = paste(
        oracle_bundle$detrendq$source_anchor,
        "+ transform_detrendq / transform_common package surface"
      ),
      seasonal_fit = oracle_bundle$detrendq$python_yhat,
      y_trans = oracle_bundle$detrendq$python_ydot,
      ydot_postavg = rep(0.25, length(oracle_bundle$detrendq$python_ydot)),
      firstpost = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
      n_pre = rep(6L, length(oracle_bundle$detrendq$python_ydot))
    )
  )
}

make_e906_forced_bundle_package_contracts <- function() {
  monthly_panel <- make_e906_monthly_panel()
  monthly_lookup <- c(
    `1` = 107,
    `2` = 102,
    `3` = 103,
    `4` = 104,
    `5` = 105,
    `6` = 106,
    `7` = 107,
    `8` = 108,
    `9` = 109,
    `10` = 110,
    `11` = 111,
    `12` = 112
  )
  monthly_seasonal_fit <- unname(monthly_lookup[as.character(monthly_panel$month)])
  monthly_y_trans <- monthly_panel$y - monthly_seasonal_fit
  monthly_firstpost <- monthly_panel$post == 1L &
    monthly_panel$tindex == min(monthly_panel$tindex[monthly_panel$post == 1L])
  monthly_n_pre <- rep(13L, nrow(monthly_panel))
  monthly_ydot_postavg <- rep(monthly_y_trans[monthly_panel$post == 1L][[1L]], nrow(monthly_panel))

  panel_data <- make_e906_panel_data()
  detrendq_y_trans <- ifelse(panel_data$D == 1L, 1.5, 0)
  detrendq_seasonal_fit <- panel_data$y - detrendq_y_trans
  detrendq_firstpost <- panel_data$post == 1L &
    panel_data$tindex == min(panel_data$tindex[panel_data$post == 1L])
  detrendq_ydot_postavg_by_id <- c(`1` = 0, `2` = 0, `3` = 1.5, `4` = 1.5)
  detrendq_ydot_postavg <- unname(detrendq_ydot_postavg_by_id[as.character(panel_data$id)])
  detrendq_n_pre <- rep(6L, nrow(panel_data))

  biannual_panel <- make_e906_biannual_identity_panel()
  biannual_firstpost <- biannual_panel$post == 1L &
    biannual_panel$tindex == min(biannual_panel$tindex[biannual_panel$post == 1L])
  biannual_postavg_by_id <- vapply(
    split(biannual_panel$expected_ydot[biannual_panel$post == 1L], biannual_panel$id[biannual_panel$post == 1L]),
    mean,
    numeric(1)
  )
  biannual_ydot_postavg <- unname(biannual_postavg_by_id[as.character(biannual_panel$id)])
  biannual_n_pre <- rep(6L, nrow(biannual_panel))

  list(
    monthly_demeanq = list(
      seasonal_fit = monthly_seasonal_fit,
      y_trans = monthly_y_trans,
      ydot_postavg = monthly_ydot_postavg,
      firstpost = monthly_firstpost,
      n_pre = monthly_n_pre
    ),
    sparse_demeanq = list(
      error = paste(
        "Unit 1 has 2 pre-period observation(s) with 2 distinct season(s).",
        "rolling('demeanq') requires at least 3 observation(s) so that df >= 1."
      )
    ),
    detrendq_alias_excluded = list(
      seasonal_fit = detrendq_seasonal_fit,
      y_trans = detrendq_y_trans,
      ydot_postavg = detrendq_ydot_postavg,
      firstpost = detrendq_firstpost,
      n_pre = detrendq_n_pre
    ),
    biannual_demeanq = list(
      seasonal_fit = biannual_panel$expected_seasonal_fit,
      y_trans = biannual_panel$expected_ydot,
      ydot_postavg = biannual_ydot_postavg,
      firstpost = biannual_firstpost,
      n_pre = biannual_n_pre
    )
  )
}

make_e906_post_summary_package_contracts <- function() {
  panel_data <- make_e906_panel_data()
  unit_effect_by_id <- c(`1` = 0.0, `2` = 0.2, `3` = 0.5, `4` = 0.7)
  unit_effect <- unname(unit_effect_by_id[as.character(panel_data$id)])
  firstpost <- panel_data$post == 1L &
    panel_data$tindex == min(panel_data$tindex[panel_data$post == 1L])

  demeanq_base_lookup <- c(`1` = 5.9, `2` = 6.6, `3` = 7.3, `4` = 8.0)
  demeanq_excluded_lookup <- c(`1` = 5.9, `2` = 6.6, `3` = 6.7, `4` = 7.4)
  demeanq_postavg_by_id <- c(`1` = 1.8, `2` = 1.8, `3` = 3.3, `4` = 3.3)

  demeanq_base_seasonal_fit <- unname(
    demeanq_base_lookup[as.character(panel_data$quarter)]
  ) + unit_effect
  demeanq_excluded_seasonal_fit <- unname(
    demeanq_excluded_lookup[as.character(panel_data$quarter)]
  ) + unit_effect
  demeanq_postavg <- unname(demeanq_postavg_by_id[as.character(panel_data$id)])

  detrendq_y_trans <- ifelse(panel_data$D == 1L, 1.5, 0)
  detrendq_postavg_by_id <- c(`1` = 0.0, `2` = 0.0, `3` = 1.5, `4` = 1.5)
  detrendq_postavg <- unname(detrendq_postavg_by_id[as.character(panel_data$id)])

  list(
    demeanq_base = list(
      source_anchor = "B6 + Task E9-06.7 transform_common(demeanq) package surface",
      seasonal_fit = demeanq_base_seasonal_fit,
      y_trans = panel_data$y - demeanq_base_seasonal_fit,
      ydot_postavg = demeanq_postavg,
      firstpost = firstpost,
      n_pre = rep(8L, nrow(panel_data))
    ),
    demeanq_excluded = list(
      source_anchor = "B6 exclude_pre_periods transform_common(demeanq) package surface",
      seasonal_fit = demeanq_excluded_seasonal_fit,
      y_trans = panel_data$y - demeanq_excluded_seasonal_fit,
      ydot_postavg = demeanq_postavg,
      firstpost = firstpost,
      n_pre = rep(6L, nrow(panel_data))
    ),
    detrendq_base = list(
      source_anchor = "C8/C9 transform_common(detrendq) package surface",
      seasonal_fit = panel_data$y - detrendq_y_trans,
      y_trans = detrendq_y_trans,
      ydot_postavg = detrendq_postavg,
      firstpost = firstpost,
      n_pre = rep(8L, nrow(panel_data))
    )
  )
}

make_e906_weekly_frequency_package_contracts <- function() {
  weekly_panel <- make_e906_weekly_detrend_panel()
  pre_panel <- weekly_panel[weekly_panel$post == 0L, , drop = FALSE]
  t_bar_pre <- mean(pre_panel$tindex)
  manual_fit <- stats::lm(
    y ~ I(tindex - t_bar_pre) + factor(week),
    data = pre_panel
  )
  weekly_seasonal_fit <- unname(stats::predict(manual_fit, newdata = weekly_panel))
  weekly_y_trans <- weekly_panel$y - weekly_seasonal_fit
  weekly_firstpost <- weekly_panel$post == 1L &
    weekly_panel$tindex == min(weekly_panel$tindex[weekly_panel$post == 1L])
  weekly_ydot_postavg <- rep(
    mean(weekly_y_trans[weekly_panel$post == 1L], na.rm = TRUE),
    nrow(weekly_panel)
  )
  frequency_detection <- lapply(
    make_e906_frequency_detection_cases(),
    function(detection_case) {
      list(
        frequency = detection_case$expected_frequency,
        Q = detection_case$expected_Q,
        confidence = 1
      )
    }
  )

  list(
    frequency_detection = frequency_detection,
    weekly_detrend = list(
      source_anchor = "C7 transform_detrendq weekly package surface",
      seasonal_fit = weekly_seasonal_fit,
      y_trans = weekly_y_trans,
      ydot_postavg = weekly_ydot_postavg,
      firstpost = weekly_firstpost,
      n_pre = rep(sum(weekly_panel$post == 0L), nrow(weekly_panel))
    ),
    invalid_detrendq_error = list(
      source_anchor = "C10 invalid seasonal dispatch surface",
      message = paste(
        "Seasonal column 'quarter' contains invalid values: 5.",
        "Expected integer values in {1, 2, ..., 4}."
      )
    )
  )
}

make_e906_validation_boundary_package_contracts <- function() {
  list(
    q_lt_two_error = list(
      source_anchor = "A5 seasonal validation exact Q boundary",
      class = "lwdid_invalid_parameter",
      message = "'Q' must be an integer >= 2. Got: 1"
    ),
    coverage_warning = list(
      source_anchor = "A2 seasonal coverage warning exact package surface",
      warnings = c(
        paste(
          "Unit 1 has post-treatment quarter(s) 3, 4 that were not observed in the",
          "pre-treatment window 1, 2. R records this as a warning and allows",
          "the seasonal gate to pass."
        )
      ),
      K = 5L,
      n_units_valid = 1L,
      n_units_invalid = 0L,
      all_covered = FALSE,
      units_with_gaps = list(c(3L, 4L)),
      total_uncovered = 2L
    ),
    exclude_pre_base = list(
      source_anchor = "A6 exclude_pre_periods baseline validation surface",
      K = 5L,
      n_units_valid = 1L,
      n_units_invalid = 0L,
      all_covered = TRUE,
      total_uncovered = 0L
    ),
    exclude_pre_shifted_error = list(
      source_anchor = "A6 exclude_pre_periods shifted validation boundary",
      class = "lwdid_insufficient_pre_periods",
      message = paste(
        "Unit 1 has 3 pre-period observation(s) with 3 distinct season(s).",
        "rolling('demeanq') requires at least 4 observation(s) so that df >= 1."
      )
    ),
    min_global_pre_error = list(
      source_anchor = "A7 min_global_pre_periods validation boundary",
      class = "lwdid_insufficient_pre_periods",
      message = "rolling('detrendq') requires at least 2 pre-treatment period(s). Found: 1."
    ),
    quarter_alias = list(
      source_anchor = "A8 quarter alias validation surface",
      K = 5L,
      n_units_valid = 1L,
      n_units_invalid = 0L,
      all_covered = TRUE,
      total_uncovered = 0L,
      warnings = character(0)
    )
  )
}

make_e906_validation_error_package_contracts <- function() {
  list(
    invalid_range_error = list(
      source_anchor = "A1 seasonal validation invalid-range error surface",
      class = "lwdid_invalid_parameter",
      message = paste(
        "Seasonal column 'quarter' contains invalid values: 5.",
        "Expected integer values in {1, 2, ..., 4}."
      )
    ),
    insufficient_pre_error = list(
      source_anchor = "A3 seasonal validation insufficient-pre error surface",
      class = "lwdid_insufficient_pre_periods",
      message = paste(
        "Unit 1 has 2 pre-period observation(s) with 2 distinct season(s).",
        "rolling('demeanq') requires at least 3 observation(s) so that df >= 1."
      )
    )
  )
}

make_e906_test_inversion_package_contracts <- function() {
  list(
    ci_surface = list(
      source_anchor = "D16 test_inversion CI exact package surface",
      ci_method = "test_inversion",
      actual_n_bootstrap = 999L,
      att = -0.054455445544589,
      ci_lower = -0.935388960667232,
      ci_upper = 0.726163419325528,
      pvalue = 0.806806806806807
    ),
    width_surface = list(
      source_anchor = "D17 test_inversion width-band exact package surface",
      test_inversion_actual_n_bootstrap = 999L,
      percentile_t_actual_n_bootstrap = 32L,
      width_ratio = 1.082421652750608,
      relative_gap = 0.082421652750608
    ),
    monotonicity_surface = list(
      source_anchor = "D19 test_inversion monotonicity exact package surface",
      att = -0.054455445544589,
      pvalue_at_att = 1,
      pvalue_far_upper = 0.03125,
      pvalue_far_lower = 0.0625
    )
  )
}

make_e906_public_boundary_package_contracts <- function() {
  list(
    public_test_inversion = list(
      source_anchor = "E7 independent public test_inversion CI exact package surface",
      ci_method = "test_inversion",
      actual_n_bootstrap = 99L,
      att = 1.21679166666667,
      pvalue = 0.131313131313131,
      ci_lower = 0.288418599727285,
      ci_upper = 2.14516473360606,
      original_se = 0.232093266734846
    ),
    controls_null_same_path = list(
      source_anchor = "F5 controls-NULL native same-path WCB exact package surface",
      actual_n_bootstrap = 256L,
      att = 0.549463375000001,
      pvalue = 0.33984375,
      ci_lower = -0.517061717711179,
      ci_upper = 1.61598846771118
    ),
    degenerate_zero_se = list(
      source_anchor = "F7 degenerate zero-SE test_inversion exact package surface",
      ci_method = "test_inversion",
      actual_n_bootstrap = 49L,
      att = -0x1.bfffffffffffdp-51,
      original_se = 0x1.ea33e2c83c13ep-51,
      pvalue = NaN,
      rejection_rate = NaN,
      ci_lower = NaN,
      ci_upper = NaN
    )
  )
}

test_that("Task E9-06.0 seasonal frequency fixtures are validation-ready for quarterly/monthly/weekly data", {
  helper_names <- c(
    "make_e906_quarterly_panel",
    "make_e906_monthly_panel",
    "make_e906_weekly_panel"
  )

  for (helper_name in helper_names) {
    expect_true(
      exists(helper_name, mode = "function"),
      info = paste("missing seasonal frequency fixture helper:", helper_name)
    )
  }

  if (!all(vapply(helper_names, exists, logical(1), mode = "function"))) {
    return(invisible(NULL))
  }

  fixture_cases <- list(
    list(
      data = make_e906_quarterly_panel(),
      season_var = "quarter",
      Q = 4L,
      expected_rows = 6L,
      expected_K = 5L
    ),
    list(
      data = make_e906_monthly_panel(),
      season_var = "month",
      Q = 12L,
      expected_rows = 14L,
      expected_K = 13L
    ),
    list(
      data = make_e906_weekly_panel(),
      season_var = "week",
      Q = 52L,
      expected_rows = 54L,
      expected_K = 53L
    )
  )

  for (fixture_case in fixture_cases) {
    expect_s3_class(fixture_case$data, "data.frame")
    expect_identical(nrow(fixture_case$data), fixture_case$expected_rows)
    expect_true(all(c("id", "tindex", "y", fixture_case$season_var, "post") %in% names(fixture_case$data)))

    validation <- lwdid:::.validate_seasonal_inputs(
      data = fixture_case$data,
      y = "y",
      season_var = fixture_case$season_var,
      Q = fixture_case$Q,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )

    expect_true(isTRUE(validation$is_valid))
    expect_identical(validation$K, fixture_case$expected_K)
    expect_true(isTRUE(validation$season_coverage$all_covered))
    expect_identical(validation$n_units_invalid, 0L)
  }
})

test_that("Task E9-06.0 deterministic seasonal fixtures preserve known-pattern metadata", {
  helper_names <- c(
    "make_e906_known_seasonal_panel",
    "make_e906_known_trend_seasonal_panel",
    "make_e906_missing_q4_panel"
  )

  for (helper_name in helper_names) {
    expect_true(
      exists(helper_name, mode = "function"),
      info = paste("missing deterministic seasonal fixture helper:", helper_name)
    )
  }

  if (!all(vapply(helper_names, exists, logical(1), mode = "function"))) {
    return(invisible(NULL))
  }

  known_seasonal <- make_e906_known_seasonal_panel()
  expect_true(all(c("id", "tindex", "quarter", "post", "y", "epsilon", "expected_seasonal_fit") %in% names(known_seasonal)))
  expect_equal(
    known_seasonal$y,
    10 + known_seasonal$expected_seasonal_fit + known_seasonal$epsilon,
    tolerance = 1e-12
  )
  expect_equal(
    known_seasonal$expected_seasonal_fit[match(1:4, known_seasonal$quarter)],
    c(0, 2, 3, 4),
    tolerance = 1e-12
  )

  known_trend <- make_e906_known_trend_seasonal_panel()
  expect_true(all(c("id", "tindex", "quarter", "post", "y", "epsilon", "expected_trend_component", "expected_seasonal_fit") %in% names(known_trend)))
  expect_equal(
    known_trend$y,
    5 + known_trend$expected_trend_component + known_trend$expected_seasonal_fit + known_trend$epsilon,
    tolerance = 1e-12
  )
  expect_equal(
    diff(known_trend$expected_trend_component),
    rep(0.3, nrow(known_trend) - 1L),
    tolerance = 1e-12
  )

  missing_q4 <- make_e906_missing_q4_panel()
  expect_false(4L %in% missing_q4$quarter[missing_q4$post == 0L])
  expect_true(4L %in% missing_q4$quarter[missing_q4$post == 1L])
})

test_that("Task E9-06.0 boundary fixtures cover single-cluster and constant-outcome edge cases", {
  helper_names <- c(
    "make_e906_single_cluster_panel",
    "make_e906_constant_outcome_panel"
  )

  for (helper_name in helper_names) {
    expect_true(
      exists(helper_name, mode = "function"),
      info = paste("missing boundary fixture helper:", helper_name)
    )
  }

  if (!all(vapply(helper_names, exists, logical(1), mode = "function"))) {
    return(invisible(NULL))
  }

  single_cluster <- make_e906_single_cluster_panel()
  expect_identical(length(unique(single_cluster$cluster)), 1L)
  expect_identical(unique(single_cluster$post), c(0L, 1L))
  expect_true(all(single_cluster$D %in% c(0, 1)))

  constant_outcome <- make_e906_constant_outcome_panel()
  expect_identical(length(unique(constant_outcome$Y)), 1L)
  expect_gt(length(unique(constant_outcome$cluster)), 1L)
  expect_true(all(constant_outcome$D %in% c(0, 1)))
})

test_that("Task E9-06.0 fixture roster covers sparse, biannual, and small-sample seasonal panels", {
  helper_names <- c(
    "make_e906_sparse_panel",
    "make_e906_constant_time_panel",
    "make_e906_all_sparse_panel",
    "make_e906_biannual_panel",
    "make_e906_small_ti_panel"
  )

  for (helper_name in helper_names) {
    expect_true(
      exists(helper_name, mode = "function"),
      info = paste("missing Task E9-06.0 fixture helper:", helper_name)
    )
  }

  if (!all(vapply(helper_names, exists, logical(1), mode = "function"))) {
    return(invisible(NULL))
  }

  sparse_panel <- make_e906_sparse_panel()
  expect_true(all(c("id", "tindex", "quarter", "post", "y") %in% names(sparse_panel)))
  expect_lt(sum(sparse_panel$post == 0L), length(unique(sparse_panel$quarter[sparse_panel$post == 0L])) + 1L)

  constant_time_panel <- make_e906_constant_time_panel()
  pre_panel <- constant_time_panel[constant_time_panel$post == 0L, , drop = FALSE]
  pre_tindex_variances <- vapply(
    split(pre_panel$tindex, pre_panel$id),
    stats::var,
    numeric(1)
  )
  expect_equal(pre_tindex_variances, c(`1` = 0, `2` = 0), tolerance = 1e-12)
  expect_identical(length(unique(pre_panel$tindex)), 2L)

  all_sparse_panel <- make_e906_all_sparse_panel()
  pre_counts_are_sparse <- vapply(
    split(all_sparse_panel, all_sparse_panel$id),
    function(unit_data) {
      pre_periods <- unit_data[unit_data$post == 0L, , drop = FALSE]
      nrow(pre_periods) <= length(unique(pre_periods$quarter))
    },
    logical(1)
  )
  expect_true(all(pre_counts_are_sparse))

  biannual_panel <- make_e906_biannual_panel()
  biannual_validation <- lwdid:::.validate_seasonal_inputs(
    data = biannual_panel,
    y = "y",
    season_var = "halfyear",
    Q = 2L,
    ivar = "id",
    tindex = "tindex",
    post = "post",
    method = "demeanq",
    exclude_pre_periods = 0L,
    min_global_pre_periods = 1L
  )
  expect_true(isTRUE(biannual_validation$is_valid))
  expect_identical(biannual_validation$K, 3L)

  small_ti_panel <- make_e906_small_ti_panel()
  expect_identical(nrow(small_ti_panel), 12L)
  expect_identical(length(unique(small_ti_panel$cluster)), 3L)
  expect_true(all(small_ti_panel$D %in% c(0, 1)))
})

test_that("Task E9-06.0 fixture roster covers panel, clustered, and collinearity helpers", {
  helper_names <- c(
    "make_e906_collinear_panel",
    "make_e906_restricted_model_panel",
    "make_e906_panel_data",
    "make_e906_seasonal_data",
    "make_e906_clustered_data",
    "make_e906_seasonal_clustered_data"
  )

  for (helper_name in helper_names) {
    expect_true(
      exists(helper_name, mode = "function"),
      info = paste("missing Task E9-06.0 integration fixture helper:", helper_name)
    )
  }

  if (!all(vapply(helper_names, exists, logical(1), mode = "function"))) {
    return(invisible(NULL))
  }

  collinear_panel <- make_e906_collinear_panel()
  expect_true(all(c("control_a", "control_b") %in% names(collinear_panel)))
  expect_equal(collinear_panel$control_b, 2 * collinear_panel$control_a, tolerance = 1e-12)

  restricted_model_panel <- make_e906_restricted_model_panel()
  expect_true(all(
    c("cluster", "control_a", "control_b", "D", "y") %in% names(restricted_model_panel)
  ))
  expect_identical(length(unique(restricted_model_panel$cluster)), 8L)
  expect_true(all(restricted_model_panel$D %in% c(0, 1)))

  panel_data <- make_e906_panel_data()
  expect_true(all(c("id", "tindex", "quarter", "post", "D", "y") %in% names(panel_data)))
  expect_gt(length(unique(panel_data$id)), 1L)

  seasonal_data <- make_e906_seasonal_data()
  expect_true(all(c("quarter", "post", "D", "y") %in% names(seasonal_data)))

  clustered_data <- make_e906_clustered_data()
  expect_true(all(c("cluster", "post", "D", "y") %in% names(clustered_data)))
  expect_gt(length(unique(clustered_data$cluster)), 1L)

  seasonal_clustered_data <- make_e906_seasonal_clustered_data()
  expect_true(all(c("quarter", "cluster", "post", "D", "y") %in% names(seasonal_clustered_data)))
  expect_gt(length(unique(seasonal_clustered_data$cluster)), 1L)
})

test_that("A5 seasonal validation rejects Q values below 2", {
  q_lt_two_panel <- data.frame(
    id = rep(1L, 4L),
    tindex = 1:4,
    quarter = rep(1L, 4L),
    post = c(0L, 0L, 1L, 1L),
    y = c(10.0, 10.5, 11.0, 11.5)
  )

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = q_lt_two_panel,
      y = "y",
      season_var = "quarter",
      Q = 1L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    ),
    regexp = "'Q' must be an integer >= 2|integer >= 2",
    class = "lwdid_invalid_parameter"
  )
})

test_that("A4 seasonal frequency auto-detection stays executable for quarterly monthly and weekly panels", {
  expect_true(
    exists("make_e906_frequency_detection_cases", mode = "function"),
    info = "missing Task E9-06.1 detection-case helper"
  )

  if (!exists("make_e906_frequency_detection_cases", mode = "function")) {
    return(invisible(NULL))
  }

  detection_cases <- make_e906_frequency_detection_cases()
  expect_equal(sort(names(detection_cases)), c("monthly", "quarterly", "weekly"))

  for (case_name in names(detection_cases)) {
    detection_case <- detection_cases[[case_name]]
    detection <- lwdid:::.auto_detect_frequency(
      data = detection_case$data,
      tvar = "year",
      ivar = "id"
    )

    expect_true(is.list(detection), info = case_name)
    expect_equal(detection$frequency, detection_case$expected_frequency, info = case_name)
    expect_equal(detection$Q, detection_case$expected_Q, info = case_name)
    expect_true(detection$confidence > 0.5, info = case_name)
  }
})

test_that("A-class seasonal validation boundary cases stay executable in the owner bundle", {
  expect_true(
    exists("make_e906_validation_case_bundle", mode = "function"),
    info = "missing Task E9-06.1 validation-case helper"
  )

  if (!exists("make_e906_validation_case_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  validation_cases <- make_e906_validation_case_bundle()

  expect_error(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a1_invalid_range$call),
    regexp = validation_cases$a1_invalid_range$regexp,
    class = "lwdid_invalid_parameter"
  )

  expect_warning(
    a2_result <- do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a2_coverage_warning$call),
    regexp = validation_cases$a2_coverage_warning$regexp,
    class = "lwdid_data"
  )
  expect_true(isTRUE(a2_result$is_valid))
  expect_false(isTRUE(a2_result$season_coverage$all_covered))
  expect_identical(unname(a2_result$season_coverage$units_with_gaps[[1L]]), c(3L, 4L))

  expect_error(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a3_insufficient_pre$call),
    regexp = validation_cases$a3_insufficient_pre$regexp,
    class = "lwdid_insufficient_pre_periods"
  )

  base_a6 <- do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a6_exclude_pre_base$call)
  expect_true(isTRUE(base_a6$is_valid))
  expect_identical(base_a6$K, 5L)
  expect_error(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a6_exclude_pre_shifted$call),
    regexp = validation_cases$a6_exclude_pre_shifted$regexp,
    class = "lwdid_insufficient_pre_periods"
  )

  expect_error(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a7_min_global_pre$call),
    regexp = validation_cases$a7_min_global_pre$regexp,
    class = "lwdid_insufficient_pre_periods"
  )

  a8_result <- do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a8_quarter_alias$call)
  expect_true(isTRUE(a8_result$is_valid))
  expect_identical(a8_result$K, 5L)
})

test_that("Task E9-06.0.3 B4 preload keeps the seasonal demeanq oracle machine-readable", {
  expect_true(
    exists("make_e906_python_reference_bundle", mode = "function"),
    info = "missing Task E9-06.0.3 Python oracle preload helper"
  )

  if (!exists("make_e906_python_reference_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  oracle_bundle <- make_e906_python_reference_bundle()
  demeanq_case <- oracle_bundle$demeanq
  demeanq_result <- lwdid:::.demeanq_unit(
    unit_data = demeanq_case$fixture,
    y = "y",
    season_var = "quarter",
    post = "post_",
    Q = 4L
  )

  expect_identical(demeanq_case$source_anchor, "tests/test_transformations.py::TestQuarterlyTransforms.test_demeanq_mve")
  expect_equal(demeanq_result$yhat, demeanq_case$python_yhat, tolerance = 1e-10)
  expect_equal(demeanq_result$ydot, demeanq_case$python_ydot, tolerance = 1e-10)
})

test_that("Task E9-06.0.3 C4 preload keeps the seasonal detrendq oracle machine-readable", {
  expect_true(
    exists("make_e906_python_reference_bundle", mode = "function"),
    info = "missing Task E9-06.0.3 Python oracle preload helper"
  )

  if (!exists("make_e906_python_reference_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  oracle_bundle <- make_e906_python_reference_bundle()
  detrendq_case <- oracle_bundle$detrendq
  detrendq_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(detrendq_case$fixture),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )

  expect_identical(detrendq_case$source_anchor, "tests/test_transformations.py::TestQuarterlyTransforms.test_detrendq_mve")
  expect_equal(detrendq_result$seasonal_fit, detrendq_case$python_yhat, tolerance = 1e-10)
  expect_equal(detrendq_result$y_trans, detrendq_case$python_ydot, tolerance = 1e-10)
})

test_that("Task E9-06.13 Layer 1 transform_demeanq keeps the Python demeanq package surface exact", {
  expect_true(
    exists("make_e906_python_transform_package_contracts", mode = "function"),
    info = "missing Task E9-06.13 helper: make_e906_python_transform_package_contracts"
  )
  expect_true(
    exists("make_e906_python_reference_bundle", mode = "function"),
    info = "missing Task E9-06.0.3 Python oracle preload helper"
  )

  if (!exists("make_e906_python_transform_package_contracts", mode = "function") ||
    !exists("make_e906_python_reference_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_python_transform_package_contracts()
  oracle_bundle <- make_e906_python_reference_bundle()
  demeanq_case <- oracle_bundle$demeanq
  demeanq_result <- lwdid:::transform_demeanq(
    dt = data.table::as.data.table(demeanq_case$fixture),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 6L,
    season_var = "quarter",
    Q = 4L,
    post = "post_"
  )

  expected <- contracts$demeanq
  expect_match(expected$source_anchor, "transform_demeanq / transform_common package surface", fixed = TRUE)
  expect_equal(demeanq_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(demeanq_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(demeanq_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(demeanq_result$firstpost, expected$firstpost)
  expect_identical(demeanq_result$n_pre, expected$n_pre)
})

test_that("Task E9-06.13 Layer 1 transform_common keeps the Python demeanq quarter alias exact", {
  expect_true(
    exists("make_e906_python_transform_package_contracts", mode = "function"),
    info = "missing Task E9-06.13 helper: make_e906_python_transform_package_contracts"
  )
  expect_true(
    exists("make_e906_python_reference_bundle", mode = "function"),
    info = "missing Task E9-06.0.3 Python oracle preload helper"
  )

  if (!exists("make_e906_python_transform_package_contracts", mode = "function") ||
    !exists("make_e906_python_reference_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_python_transform_package_contracts()
  oracle_bundle <- make_e906_python_reference_bundle()
  demeanq_case <- oracle_bundle$demeanq
  direct_result <- lwdid:::transform_demeanq(
    dt = data.table::as.data.table(demeanq_case$fixture),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 6L,
    season_var = "quarter",
    Q = 4L,
    post = "post_"
  )
  alias_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(demeanq_case$fixture),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 6L,
    rolling = "demeanq",
    quarter = "quarter",
    Q = 4L,
    post = "post_"
  )

  expected <- contracts$demeanq
  expect_equal(alias_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(alias_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(alias_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(alias_result$firstpost, expected$firstpost)
  expect_identical(alias_result$n_pre, expected$n_pre)
  expect_equal(alias_result$seasonal_fit, direct_result$seasonal_fit, tolerance = 1e-10)
  expect_equal(alias_result$y_trans, direct_result$y_trans, tolerance = 1e-10)
  expect_identical(alias_result$firstpost, direct_result$firstpost)
})

test_that("Task E9-06.13 Layer 1 transform_common keeps the Python detrendq quarter alias exact", {
  expect_true(
    exists("make_e906_python_transform_package_contracts", mode = "function"),
    info = "missing Task E9-06.13 helper: make_e906_python_transform_package_contracts"
  )
  expect_true(
    exists("make_e906_python_reference_bundle", mode = "function"),
    info = "missing Task E9-06.0.3 Python oracle preload helper"
  )

  if (!exists("make_e906_python_transform_package_contracts", mode = "function") ||
    !exists("make_e906_python_reference_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_python_transform_package_contracts()
  oracle_bundle <- make_e906_python_reference_bundle()
  detrendq_case <- oracle_bundle$detrendq
  direct_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(detrendq_case$fixture),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )
  alias_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(detrendq_case$fixture),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    rolling = "detrendq",
    quarter = "quarter",
    Q = 4L,
    post = "post"
  )

  expected <- contracts$detrendq
  expect_match(expected$source_anchor, "transform_detrendq / transform_common package surface", fixed = TRUE)
  expect_equal(alias_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(alias_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(alias_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(alias_result$firstpost, expected$firstpost)
  expect_identical(alias_result$n_pre, expected$n_pre)
  expect_equal(alias_result$seasonal_fit, direct_result$seasonal_fit, tolerance = 1e-10)
  expect_equal(alias_result$y_trans, direct_result$y_trans, tolerance = 1e-10)
  expect_identical(alias_result$firstpost, direct_result$firstpost)
})

test_that("B5 monthly demeanq keeps every transformed period finite and centered on the pre window", {
  monthly_panel <- make_e906_monthly_panel()

  demeanq_result <- lwdid:::.demeanq_unit(
    unit_data = monthly_panel,
    y = "y",
    season_var = "month",
    post = "post",
    Q = 12L
  )

  expect_false(any(is.na(demeanq_result$ydot)))
  expect_lt(abs(mean(demeanq_result$ydot[monthly_panel$post == 0L], na.rm = TRUE)), 1e-10)
  expect_equal(unname(demeanq_result$ydot[monthly_panel$post == 1L]), 12)
})

test_that("B1 known seasonal fixtures recover the seasonal fit and residual epsilon", {
  known_panel <- make_e906_known_seasonal_panel()

  demeanq_result <- lwdid:::.demeanq_unit(
    unit_data = known_panel,
    y = "y",
    season_var = "quarter",
    post = "post",
    Q = 4L
  )

  expect_false(any(is.na(demeanq_result$ydot)))
  expect_equal(
    demeanq_result$yhat,
    10 + known_panel$expected_seasonal_fit,
    tolerance = 1e-8
  )
  expect_equal(demeanq_result$ydot, known_panel$epsilon, tolerance = 1e-8)
})

test_that("B2 unseen post seasons fall back to the intercept-only pre-season fit", {
  missing_q4_panel <- make_e906_missing_q4_panel()
  q4_post_idx <- missing_q4_panel$quarter == 4L & missing_q4_panel$post == 1L
  q1_pre_mean <- mean(
    missing_q4_panel$y[missing_q4_panel$post == 0L & missing_q4_panel$quarter == 1L]
  )

  demeanq_result <- suppressWarnings(
    lwdid:::.demeanq_unit(
      unit_data = missing_q4_panel,
      y = "y",
      season_var = "quarter",
      post = "post",
      Q = 4L
    )
  )

  expect_true(any(q4_post_idx))
  expect_false(any(is.na(demeanq_result$ydot[q4_post_idx])))
  expect_equal(
    demeanq_result$yhat[q4_post_idx],
    rep(q1_pre_mean, sum(q4_post_idx)),
    tolerance = 1e-10
  )
})

test_that("B3 sparse seasonal units warn and return all-NA demeanq outputs", {
  sparse_panel <- make_e906_sparse_panel()
  demeanq_result <- NULL

  expect_warning(
    demeanq_result <- lwdid:::.demeanq_unit(
      unit_data = sparse_panel,
      y = "y",
      season_var = "quarter",
      post = "post",
      Q = 4L
    ),
    "valid pre-treatment observation\\(s\\)|seasonal transformed outcomes are set to NA"
  )

  expect_true(all(is.na(demeanq_result$yhat)))
  expect_true(all(is.na(demeanq_result$ydot)))
  expect_identical(demeanq_result$n_pre, 2L)
})

test_that("B6 exclude_pre_periods changes seasonal demeanq fits and transformed outcomes", {
  panel_data <- make_e906_panel_data()

  base_result <- lwdid:::.demeanq_transform(
    data = panel_data,
    y = "y",
    ivar = "id",
    tindex = "tindex",
    season_var = "quarter",
    post = "post",
    Q = 4L,
    exclude_pre_periods = 0L
  )
  shifted_result <- lwdid:::.demeanq_transform(
    data = panel_data,
    y = "y",
    ivar = "id",
    tindex = "tindex",
    season_var = "quarter",
    post = "post",
    Q = 4L,
    exclude_pre_periods = 2L
  )

  expect_true(any(base_result$n_pre > shifted_result$n_pre))
  expect_false(isTRUE(all.equal(base_result$seasonal_fit, shifted_result$seasonal_fit)))
  expect_false(isTRUE(all.equal(base_result$y_trans, shifted_result$y_trans)))
})

test_that("Task E9-06.7 demeanq broadcasts ydot_postavg as a unit-level post mean", {
  panel_data <- make_e906_panel_data()

  demeanq_result <- lwdid:::.demeanq_transform(
    data = panel_data,
    y = "y",
    ivar = "id",
    tindex = "tindex",
    season_var = "quarter",
    post = "post",
    Q = 4L
  )

  expect_true(all(c("ydot_postavg", "firstpost") %in% names(demeanq_result)))

  for (uid in unique(panel_data$id)) {
    unit_mask <- panel_data$id == uid
    unit_postavg <- unique(demeanq_result$ydot_postavg[unit_mask])
    post_ydot <- demeanq_result$y_trans[unit_mask & panel_data$post == 1L]

    expect_identical(length(unit_postavg), 1L)
    expect_equal(unit_postavg, mean(post_ydot, na.rm = TRUE), tolerance = 1e-10)
  }
})

test_that("Task E9-06.7 demeanq firstpost marks the first original post period once per unit", {
  panel_data <- make_e906_panel_data()
  tpost1 <- min(panel_data$tindex[panel_data$post == 1L])

  demeanq_result <- lwdid:::.demeanq_transform(
    data = panel_data,
    y = "y",
    ivar = "id",
    tindex = "tindex",
    season_var = "quarter",
    post = "post",
    Q = 4L
  )

  expect_true(all(c("ydot_postavg", "firstpost") %in% names(demeanq_result)))

  firstpost_rows <- demeanq_result[demeanq_result$firstpost, , drop = FALSE]
  expect_identical(nrow(firstpost_rows), length(unique(panel_data$id)))
  expect_true(all(firstpost_rows$post == 1L))
  expect_true(all(!is.na(firstpost_rows$ydot_postavg)))

  for (uid in unique(panel_data$id)) {
    unit_firstpost <- firstpost_rows[firstpost_rows$id == uid, , drop = FALSE]
    expect_identical(nrow(unit_firstpost), 1L)
    expect_identical(unit_firstpost$tindex[[1L]], tpost1)
  }
})

test_that("C1 known trend-seasonal panel matches a centered seasonal OLS detrend fit", {
  known_panel <- make_e906_known_trend_seasonal_panel()
  pre_panel <- known_panel[known_panel$post == 0L, , drop = FALSE]
  t_bar_pre <- mean(pre_panel$tindex)
  manual_fit <- stats::lm(
    y ~ I(tindex - t_bar_pre) + factor(quarter),
    data = pre_panel
  )
  manual_yhat <- unname(stats::predict(manual_fit, newdata = known_panel))

  detrendq_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(known_panel),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 7L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )

  expect_false(any(is.na(detrendq_result$y_trans)))
  expect_equal(detrendq_result$seasonal_fit, manual_yhat, tolerance = 1e-10)
  expect_equal(detrendq_result$y_trans, known_panel$y - manual_yhat, tolerance = 1e-10)
})

test_that("C2 detrendq fit is invariant to translating the time index with g", {
  known_panel <- make_e906_known_trend_seasonal_panel()
  shifted_panel <- known_panel
  shifted_panel$tindex <- shifted_panel$tindex + 100L

  base_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(known_panel),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 7L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )
  shifted_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(shifted_panel),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 107L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )

  expect_equal(shifted_result$seasonal_fit, base_result$seasonal_fit, tolerance = 1e-10)
  expect_equal(shifted_result$y_trans, base_result$y_trans, tolerance = 1e-10)
  expect_equal(shifted_result$slope, base_result$slope, tolerance = 1e-10)
  expect_equal(
    shifted_result$t_bar_pre - base_result$t_bar_pre,
    rep(100, nrow(base_result)),
    tolerance = 1e-10
  )
})

test_that("C3 detrendq warns and returns NA outputs when each unit has zero pre-period time variance", {
  constant_time_panel <- make_e906_constant_time_panel()

  expect_warning(
    detrendq_result <- lwdid:::transform_detrendq(
      dt = data.table::as.data.table(constant_time_panel),
      y = "y",
      ivar = "id",
      tvar = "tindex",
      g = 6L,
      season_var = "quarter",
      Q = 4L,
      post = "post"
    ),
    "degenerate pre-treatment time variation"
  )

  expect_true(all(is.na(detrendq_result$seasonal_fit)))
  expect_true(all(is.na(detrendq_result$y_trans)))
  expect_true(all(is.na(detrendq_result$intercept_c)))
  expect_true(all(is.na(detrendq_result$slope)))
  expect_true(all(is.na(detrendq_result$t_bar_pre)))
})

test_that("C5 transform_common dispatches detrendq through the same seasonal-transform path", {
  panel_data <- make_e906_panel_data()

  direct_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )
  dispatch_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )

  expect_equal(dispatch_result$seasonal_fit, direct_result$seasonal_fit, tolerance = 1e-10)
  expect_equal(dispatch_result$y_trans, direct_result$y_trans, tolerance = 1e-10)
  expect_equal(dispatch_result$ydot_postavg, direct_result$ydot_postavg, tolerance = 1e-10)
  expect_identical(dispatch_result$firstpost, direct_result$firstpost)
})

test_that("C6 detrendq dispatch errors with an example when seasonal metadata is missing", {
  panel_data <- make_e906_panel_data()

  expect_error(
    lwdid:::transform_common(
      dt = data.table::as.data.table(panel_data),
      y = "y",
      ivar = "id",
      tvar = "tindex",
      g = 9L,
      rolling = "detrendq",
      post = "post",
      Q = 4L
    ),
    "rolling='detrendq'.*season_var='quarter'.*Q=4"
  )
})

test_that("C7 weekly detrendq stays finite once the pre window covers 52 seasonal buckets plus slope df", {
  weekly_panel <- make_e906_weekly_detrend_panel()
  pre_panel <- weekly_panel[weekly_panel$post == 0L, , drop = FALSE]
  t_bar_pre <- mean(pre_panel$tindex)
  manual_fit <- stats::lm(
    y ~ I(tindex - t_bar_pre) + factor(week),
    data = pre_panel
  )
  manual_yhat <- unname(stats::predict(manual_fit, newdata = weekly_panel))

  detrendq_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(weekly_panel),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 55L,
    season_var = "week",
    Q = 52L,
    post = "post"
  )

  expect_false(any(is.na(detrendq_result$y_trans)))
  expect_equal(detrendq_result$seasonal_fit, manual_yhat, tolerance = 1e-10)
  expect_equal(detrendq_result$y_trans, weekly_panel$y - manual_yhat, tolerance = 1e-10)
})

test_that("C8 detrendq broadcasts ydot_postavg as a unit-level post mean", {
  panel_data <- make_e906_panel_data()

  detrendq_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )

  expect_true(all(c("ydot_postavg", "firstpost") %in% names(detrendq_result)))

  for (uid in unique(panel_data$id)) {
    unit_mask <- panel_data$id == uid
    unit_postavg <- unique(detrendq_result$ydot_postavg[unit_mask])
    post_ydot <- detrendq_result$y_trans[unit_mask & panel_data$post == 1L]

    expect_identical(length(unit_postavg), 1L)
    expect_equal(unit_postavg, mean(post_ydot, na.rm = TRUE), tolerance = 1e-10)
  }
})

test_that("C9 detrendq firstpost marks the first original post period once per unit", {
  panel_data <- make_e906_panel_data()
  tpost1 <- min(panel_data$tindex[panel_data$post == 1L])

  detrendq_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    season_var = "quarter",
    Q = 4L,
    post = "post"
  )

  expect_true(all(c("ydot_postavg", "firstpost") %in% names(detrendq_result)))

  firstpost_rows <- detrendq_result[detrendq_result$firstpost, , drop = FALSE]
  expect_identical(nrow(firstpost_rows), length(unique(panel_data$id)))
  expect_true(all(firstpost_rows$post == 1L))
  expect_true(all(!is.na(firstpost_rows$ydot_postavg)))

  for (uid in unique(panel_data$id)) {
    unit_firstpost <- firstpost_rows[firstpost_rows$id == uid, , drop = FALSE]
    expect_identical(nrow(unit_firstpost), 1L)
    expect_identical(unit_firstpost$tindex[[1L]], tpost1)
  }
})

test_that("C11 detrendq keeps ydot_postavg anchored to the original post window under exclude_pre_periods", {
  panel_data <- make_e906_panel_data()
  effective_post <- lwdid:::.compute_effective_seasonal_post(
    data = panel_data,
    ivar = "id",
    tindex = "tindex",
    post = "post",
    exclude_pre_periods = 2L
  )

  detrendq_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    season_var = "quarter",
    Q = 4L,
    post = "post",
    exclude_pre_periods = 2L
  )

  expect_true(all(c("ydot_postavg", "firstpost") %in% names(detrendq_result)))
  expect_true(any(panel_data$post == 0L & effective_post == 1L))

  unit_ids <- unique(panel_data$id)
  observed_postavg <- vapply(
    unit_ids,
    function(uid) unique(detrendq_result$ydot_postavg[panel_data$id == uid]),
    numeric(1)
  )
  original_postavg <- vapply(
    unit_ids,
    function(uid) {
      mean(
        detrendq_result$y_trans[panel_data$id == uid & panel_data$post == 1L],
        na.rm = TRUE
      )
    },
    numeric(1)
  )
  effective_postavg <- vapply(
    unit_ids,
    function(uid) {
      mean(
        detrendq_result$y_trans[panel_data$id == uid & effective_post == 1L],
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  expect_true(any(abs(original_postavg - effective_postavg) > 1e-10))
  expect_equal(observed_postavg, original_postavg, tolerance = 1e-10)
})

test_that("C10 detrendq dispatch rejects invalid seasonal values before fitting", {
  panel_data <- make_e906_panel_data()
  panel_data$quarter[[1L]] <- 5L

  expect_error(
    lwdid:::transform_common(
      dt = data.table::as.data.table(panel_data),
      y = "y",
      ivar = "id",
      tvar = "tindex",
      g = 9L,
      rolling = "detrendq",
      post = "post",
      season_var = "quarter",
      Q = 4L
    ),
    "contains invalid values: 5"
  )
})

test_that("C12 detrendq dispatch keeps the quarter alias aligned with season_var", {
  panel_data <- make_e906_panel_data()

  season_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )
  quarter_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    rolling = "detrendq",
    post = "post",
    quarter = "quarter",
    Q = 4L
  )

  expect_equal(quarter_result$seasonal_fit, season_result$seasonal_fit, tolerance = 1e-10)
  expect_equal(quarter_result$y_trans, season_result$y_trans, tolerance = 1e-10)
  expect_equal(quarter_result$ydot_postavg, season_result$ydot_postavg, tolerance = 1e-10)
  expect_identical(quarter_result$firstpost, season_result$firstpost)
})

test_that("F3 biannual demeanq keeps the minimal-Q seasonal surface finite", {
  expect_true(
    exists("make_e906_biannual_identity_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_biannual_identity_panel"
  )

  if (!exists("make_e906_biannual_identity_panel", mode = "function")) {
    return(invisible(NULL))
  }

  biannual_panel <- make_e906_biannual_identity_panel()
  biannual_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(biannual_panel),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 7L,
    rolling = "demeanq",
    post = "post",
    season_var = "halfyear",
    Q = 2L
  )

  expect_true(all(is.finite(biannual_result$y_trans)))
  expect_equal(
    biannual_result$seasonal_fit,
    biannual_panel$expected_seasonal_fit,
    tolerance = 1e-10
  )
  expect_equal(
    biannual_result$y_trans,
    biannual_panel$expected_ydot,
    tolerance = 1e-10
  )
  expect_identical(unique(biannual_result$n_pre), 6L)
})

test_that("F1 single-cluster native WCB returns NaN inference surfaces", {
  single_cluster_panel <- make_e906_single_cluster_panel()
  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = single_cluster_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(wcb_result, "lwdid_wcb_result")
  expect_true(is.nan(wcb_result$original_se))
  expect_true(is.nan(wcb_result$pvalue))
  expect_true(is.nan(wcb_result$rejection_rate))
  expect_true(is.nan(wcb_result$ci_lower))
  expect_true(is.nan(wcb_result$ci_upper))
})

test_that("broader hardening native WCB treats constant-treatment designs as degenerate without controls", {
  constant_treatment_panel <- make_e906_constant_treatment_cluster_panel()

  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = constant_treatment_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(wcb_result, "lwdid_wcb_result")
  expect_true(is.nan(wcb_result$original_se))
  expect_true(is.nan(wcb_result$pvalue))
  expect_true(is.nan(wcb_result$rejection_rate))
  expect_true(is.nan(wcb_result$ci_lower))
  expect_true(is.nan(wcb_result$ci_upper))
  expect_true(is.nan(wcb_result$t_stat_original))
  expect_identical(wcb_result$n_clusters, 4L)
})

test_that("broader hardening native WCB treats constant-treatment designs as degenerate with controls", {
  constant_treatment_panel <- make_e906_constant_treatment_cluster_panel()

  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = constant_treatment_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = c("control_a", "control_b"),
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      restricted_model = "with_controls",
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(wcb_result, "lwdid_wcb_result")
  expect_true(is.nan(wcb_result$original_se))
  expect_true(is.nan(wcb_result$pvalue))
  expect_true(is.nan(wcb_result$rejection_rate))
  expect_true(is.nan(wcb_result$ci_lower))
  expect_true(is.nan(wcb_result$ci_upper))
  expect_true(is.nan(wcb_result$t_stat_original))
  expect_identical(wcb_result$n_clusters, 4L)
})

test_that("F2 all-invalid seasonal demean returns NaN transformed outcomes at the package surface", {
  all_sparse_panel <- make_e906_all_sparse_panel()
  demeanq_result <- NULL

  expect_warning(
    demeanq_result <- lwdid:::.demeanq_transform(
      data = all_sparse_panel,
      y = "y",
      ivar = "id",
      tindex = "tindex",
      season_var = "quarter",
      post = "post",
      Q = 4L
    ),
    "valid pre-treatment observation\\(s\\)|seasonal transformed outcomes are set to NA"
  )

  expect_true(all(is.nan(demeanq_result$seasonal_fit)))
  expect_true(all(is.nan(demeanq_result$y_trans)))
  expect_true(all(is.nan(demeanq_result$ydot_postavg)))
  expect_false(any(demeanq_result$firstpost))
})

test_that("F4 G=12 Rademacher full enumeration expands to 4096 sign patterns", {
  weight_info <- lwdid:::.resolve_wild_bootstrap_weights(
    n_clusters = 12L,
    requested_n_bootstrap = 999L,
    weight_type = "rademacher",
    full_enumeration = NULL
  )

  expect_true(isTRUE(weight_info$full_enumeration))
  expect_identical(dim(weight_info$weights), c(4096L, 12L))
  expect_identical(weight_info$actual_n_bootstrap, 4096L)
  expect_equal(sort(unique(as.vector(weight_info$weights))), c(-1, 1))
})

test_that("F4 public native WCB keeps the G=12 boundary on the 4096-draw exact path", {
  expect_true(
    exists("make_e906_public_g12_full_enumeration_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_g12_full_enumeration_panel"
  )

  if (!exists("make_e906_public_g12_full_enumeration_panel", mode = "function")) {
    return(invisible(NULL))
  }

  public_g12_panel <- make_e906_public_g12_full_enumeration_panel()
  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = public_g12_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 199L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_true(isTRUE(wcb_result$full_enumeration))
  expect_identical(wcb_result$requested_n_bootstrap, 199L)
  expect_identical(wcb_result$actual_n_bootstrap, 4096L)
  expect_length(wcb_result$att_bootstrap, 4096L)
  expect_length(wcb_result$t_stats_bootstrap, 4096L)
})

test_that("F8 demeanq computes transformed outcomes for all periods in the biannual identity panel", {
  expect_true(
    exists("make_e906_biannual_identity_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_biannual_identity_panel"
  )

  if (!exists("make_e906_biannual_identity_panel", mode = "function")) {
    return(invisible(NULL))
  }

  biannual_panel <- make_e906_biannual_identity_panel()
  biannual_result <- lwdid:::.demeanq_transform(
    data = biannual_panel,
    y = "y",
    ivar = "id",
    tindex = "tindex",
    season_var = "halfyear",
    post = "post",
    Q = 2L
  )

  expect_false(any(is.na(biannual_result$y_trans)))
  expect_equal(
    biannual_result$y_trans,
    biannual_panel$expected_ydot,
    tolerance = 1e-10
  )
  expect_true(all(is.finite(biannual_result$y_trans[biannual_panel$post == 0L])))
  expect_true(all(is.finite(biannual_result$y_trans[biannual_panel$post == 1L])))
})

test_that("F9 demeanq pre-period residual means stay numerically zero by unit", {
  expect_true(
    exists("make_e906_biannual_identity_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_biannual_identity_panel"
  )

  if (!exists("make_e906_biannual_identity_panel", mode = "function")) {
    return(invisible(NULL))
  }

  biannual_panel <- make_e906_biannual_identity_panel()
  biannual_result <- lwdid:::.demeanq_transform(
    data = biannual_panel,
    y = "y",
    ivar = "id",
    tindex = "tindex",
    season_var = "halfyear",
    post = "post",
    Q = 2L
  )

  pre_means <- vapply(
    split(
      biannual_result$y_trans[biannual_panel$post == 0L],
      biannual_panel$id[biannual_panel$post == 0L]
    ),
    mean,
    numeric(1)
  )

  expect_lt(max(abs(pre_means)), 1e-10)
  expect_equal(unname(pre_means), c(0, 0), tolerance = 1e-10)
})

test_that("D10 Python-oracle preload keeps the intercept-only full-enumeration comparator executable", {
  expect_true(
    exists("make_e906_python_reference_bundle", mode = "function"),
    info = "missing Task E9-06.0.3 Python oracle preload helper"
  )

  if (!exists("make_e906_python_reference_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  oracle_bundle <- make_e906_python_reference_bundle()
  wcb_case <- oracle_bundle$wild_bootstrap
  wcb_result <- suppressWarnings(
    suppressMessages(
      lwdid::wild_cluster_bootstrap(
        data = wcb_case$fixture,
        y_transformed = "y_transformed",
        d = "d",
        cluster_var = "cluster",
        controls = NULL,
        n_bootstrap = 199L,
        weight_type = "rademacher",
        impose_null = TRUE,
        alpha = 0.05,
        seed = 42L,
        full_enumeration = TRUE,
        use_fwildclusterboot = FALSE,
        restricted_model = "intercept_only"
      )
    )
  )

  expect_identical(wcb_case$source_anchor, "tests/inference/test_wild_bootstrap_optimization.py::_run_bootstrap")
  expect_identical(wcb_result$actual_n_bootstrap, wcb_case$python_n_bootstrap)
  expect_identical(wcb_result$n_clusters, wcb_case$n_clusters)
  expect_true(isTRUE(wcb_result$full_enumeration))
  expect_equal(wcb_result$att, wcb_case$python_att, tolerance = 1e-10)
  expect_equal(wcb_result$t_stat_original, wcb_case$python_t_stat_original, tolerance = 1e-10)
  expect_equal(wcb_result$ci_lower, wcb_case$python_ci_lower, tolerance = 1e-10)
  expect_equal(wcb_result$ci_upper, wcb_case$python_ci_upper, tolerance = 1e-10)
  expect_equal(wcb_result$pvalue, wcb_case$python_pvalue, tolerance = 4 / wcb_case$python_n_bootstrap)
})

test_that("D10b Python multiscenario preload keeps ATT and t-stat invariants executable", {
  expect_true(
    exists("make_e906_python_reference_bundle", mode = "function"),
    info = "missing Task E9-06.0.3 Python oracle preload helper"
  )

  if (!exists("make_e906_python_reference_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  oracle_bundle <- make_e906_python_reference_bundle()
  multiscenario_case <- oracle_bundle$wild_bootstrap_multiscenario

  expect_type(multiscenario_case, "list")
  expect_length(multiscenario_case$scenarios, 6L)

  scenario_names <- names(multiscenario_case$scenarios)
  expect_equal(
    sort(scenario_names),
    sort(c(
      "impose_null_true_rademacher",
      "impose_null_true_mammen",
      "impose_null_true_webb",
      "impose_null_false_rademacher",
      "impose_null_false_mammen",
      "impose_null_false_webb"
    ))
  )

  for (scenario_name in scenario_names) {
    scenario_reference <- multiscenario_case$scenarios[[scenario_name]]

    wcb_result <- suppressWarnings(
      suppressMessages(
        lwdid::wild_cluster_bootstrap(
          data = multiscenario_case$fixture,
          y_transformed = "y_transformed",
          d = "d",
          cluster_var = "cluster",
          controls = NULL,
          n_bootstrap = multiscenario_case$n_bootstrap,
          weight_type = scenario_reference$weight_type,
          impose_null = scenario_reference$impose_null,
          alpha = multiscenario_case$alpha,
          seed = multiscenario_case$seed,
          full_enumeration = FALSE,
          use_fwildclusterboot = FALSE,
          restricted_model = "intercept_only"
        )
      )
    )

    expect_identical(wcb_result$requested_n_bootstrap, multiscenario_case$n_bootstrap)
    expect_identical(wcb_result$actual_n_bootstrap, multiscenario_case$n_bootstrap)
    expect_false(isTRUE(wcb_result$full_enumeration))
    expect_identical(wcb_result$n_clusters, multiscenario_case$n_clusters)
    expect_equal(wcb_result$att, multiscenario_case$shared_att, tolerance = 1e-10)
    expect_equal(wcb_result$att, scenario_reference$att, tolerance = 1e-10)
    expect_equal(wcb_result$t_stat_original, multiscenario_case$shared_t_stat_original, tolerance = 1e-10)
    expect_equal(wcb_result$t_stat_original, scenario_reference$t_stat_original, tolerance = 1e-10)
  }
})

test_that("D10b Python multiscenario preload keeps scenario-specific sampled bootstrap outputs finite and non-collapsed", {
  expect_true(
    exists("make_e906_python_reference_bundle", mode = "function"),
    info = "missing Task E9-06.0.3 Python oracle preload helper"
  )

  if (!exists("make_e906_python_reference_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  oracle_bundle <- make_e906_python_reference_bundle()
  multiscenario_case <- oracle_bundle$wild_bootstrap_multiscenario
  expect_type(multiscenario_case, "list")
  expect_identical(
    multiscenario_case$source_anchor,
    "tests/inference/test_wild_bootstrap_optimization.py::SCENARIOS + reference_data/wild_bootstrap_reference.json"
  )
  expect_identical(multiscenario_case$n_bootstrap, 999L)
  expect_identical(multiscenario_case$seed, 42L)
  expect_identical(multiscenario_case$n_clusters, 10L)

  scenario_reference <- multiscenario_case$scenarios
  sampled_surface <- lapply(names(scenario_reference), function(scenario_name) {
    reference <- scenario_reference[[scenario_name]]
    wcb_result <- suppressWarnings(
      suppressMessages(
        lwdid::wild_cluster_bootstrap(
          data = multiscenario_case$fixture,
          y_transformed = "y_transformed",
          d = "d",
          cluster_var = "cluster",
          controls = NULL,
          n_bootstrap = multiscenario_case$n_bootstrap,
          weight_type = reference$weight_type,
          impose_null = reference$impose_null,
          alpha = multiscenario_case$alpha,
          seed = multiscenario_case$seed,
          full_enumeration = FALSE,
          use_fwildclusterboot = FALSE,
          restricted_model = "intercept_only"
        )
      )
    )

    expect_identical(wcb_result$requested_n_bootstrap, multiscenario_case$n_bootstrap)
    expect_identical(wcb_result$actual_n_bootstrap, multiscenario_case$n_bootstrap)
    expect_false(isTRUE(wcb_result$full_enumeration))
    expect_identical(wcb_result$n_clusters, multiscenario_case$n_clusters)
    expect_true(is.finite(wcb_result$se_bootstrap), info = scenario_name)
    expect_true(is.finite(wcb_result$pvalue), info = scenario_name)
    expect_true(is.finite(wcb_result$ci_lower), info = scenario_name)
    expect_true(is.finite(wcb_result$ci_upper), info = scenario_name)
    expect_true(wcb_result$se_bootstrap > 0, info = scenario_name)
    expect_true(wcb_result$pvalue >= 0, info = scenario_name)
    expect_true(wcb_result$pvalue <= 1, info = scenario_name)
    expect_true(wcb_result$ci_lower < wcb_result$ci_upper, info = scenario_name)

    c(
      se_bootstrap = unname(wcb_result$se_bootstrap),
      pvalue = unname(wcb_result$pvalue),
      ci_lower = unname(wcb_result$ci_lower),
      ci_upper = unname(wcb_result$ci_upper)
    )
  })

  sampled_surface <- do.call(rbind, sampled_surface)
  rownames(sampled_surface) <- names(scenario_reference)

  for (scenario_name in names(scenario_reference)) {
    reference <- scenario_reference[[scenario_name]]
    expect_true(is.finite(reference$se_bootstrap), info = scenario_name)
    expect_true(is.finite(reference$pvalue), info = scenario_name)
    expect_true(is.finite(reference$ci_lower), info = scenario_name)
    expect_true(is.finite(reference$ci_upper), info = scenario_name)
  }

  expect_gt(length(unique(signif(sampled_surface[, "se_bootstrap"], 8))), 1L)
  expect_gt(length(unique(signif(sampled_surface[, "pvalue"], 8))), 1L)
  expect_gt(length(unique(signif(sampled_surface[, "ci_lower"], 8))), 1L)
  expect_gt(length(unique(signif(sampled_surface[, "ci_upper"], 8))), 1L)
})

test_that("Task E9-06.0 clustered DGP helpers materialize package-ready zero/strong-effect fixtures", {
  zero_effect <- simulate_zero_effect_clustered(N = 200L, G = 10L, seed = 42L)
  strong_effect <- simulate_strong_effect_clustered(N = 200L, G = 10L, tau = 2, seed = 42L)

  for (fixture in list(zero_effect, strong_effect)) {
    expect_s3_class(fixture, "data.frame")
    expect_identical(names(fixture), c("Y", "D", "cluster"))
    expect_identical(nrow(fixture), 200L)
    expect_identical(length(unique(fixture$cluster)), 10L)
    expect_true(all(fixture$D %in% c(0, 1)))
    expect_true(any(fixture$D == 0))
    expect_true(any(fixture$D == 1))
  }

  expect_gt(mean(strong_effect$Y[strong_effect$D == 1]), mean(zero_effect$Y[zero_effect$D == 1]))
})

test_that("D11 size frontier stays executable in the seasonal-WCB package slice", {
  scenario <- make_seasonal_wcb_decision_scenario(
    scenario_id = "story_e9_06_null_effect_non_rejection",
    true_tau = 0,
    n_simulations = 100L,
    metric_name = "non_rejection_rate",
    decision_direction = "gt"
  )

  result <- lwdid:::.run_clustering_monte_carlo_wild_bootstrap_decision_scenario(
    scenario = scenario,
    shared_dgp = make_seasonal_wcb_shared_dgp()
  )

  expect_identical(result$metric_name, "non_rejection_rate")
  expect_identical(result$n_simulations, 100L)
  expect_gte(result$metric_hit_count, 85L)
  expect_gte(result$metric_value, 0.85)
  expect_lte(result$metric_value, 1)
  expect_true(is.finite(result$metric_confidence_interval$lower))
  expect_true(is.finite(result$metric_confidence_interval$upper))
  expect_lte(result$metric_confidence_interval$lower, result$metric_value)
  expect_gte(result$metric_confidence_interval$upper, result$metric_value)
  expect_true(isTRUE(result$full_enumeration))
  expect_identical(result$actual_n_bootstrap, 1024L)
})

test_that("D12 power frontier stays executable in the seasonal-WCB package slice", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics()

  expect_identical(diagnostics$case_id, "TC-9.4.19")
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
  expect_equal(diagnostics$threshold, 0.70, tolerance = 1e-12)
  expect_equal(diagnostics$required_power_to_clear, 0.72, tolerance = 1e-12)
  expect_identical(diagnostics$result$rejection_count, 16L)
  expect_equal(diagnostics$result$metric_value, 0.32, tolerance = 1e-12)
  expect_identical(diagnostics$rejection_count_shortfall, 20L)
  expect_equal(diagnostics$threshold_gap, -0.38, tolerance = 1e-12)
  expect_false(isTRUE(diagnostics$threshold_passed))
})

test_that("D12 case-combo mechanism bridge stays executable in the seasonal-WCB package slice", {
  bridge <- lwdid:::.story_local_wcb_case_combo_mechanism_bridge()

  expect_identical(bridge$case_id, "TC-9.4.19")
  expect_identical(
    bridge$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(bridge$top_bucket_order, c(5L, 6L, 4L))
  expect_identical(bridge$low_signal_plus_elevated_se_buckets, c(5L, 6L))
  expect_identical(bridge$low_signal_dominant_buckets, 4L)
  expect_identical(bridge$residual_bucket_support, c(2L, 3L, 8L, 9L))
  expect_equal(bridge$summary$top_bucket_case_share, 0.75, tolerance = 1e-12)
  expect_equal(
    bridge$summary$low_signal_plus_elevated_se_case_share,
    0.65,
    tolerance = 1e-12
  )
  expect_equal(
    bridge$summary$low_signal_dominant_case_share,
    0.10,
    tolerance = 1e-12
  )
  expect_equal(bridge$summary$residual_case_share, 0.25, tolerance = 1e-12)
  expect_true(isTRUE(
    bridge$summary$top_bucket_family_covers_threshold_combo_majority
  ))
})

test_that("D12 threshold-crossing partition stays executable in the seasonal-WCB package slice", {
  partition <- lwdid:::.story_local_wcb_case_combo_residual_partition()

  expect_identical(partition$case_id, "TC-9.4.19")
  expect_identical(
    partition$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(partition$minimum_case_rank_to_clear_threshold, 20L)
  expect_identical(partition$overlap_bucket_support, c(2L, 3L, 4L, 5L, 6L))
  expect_identical(partition$minimum_only_bucket_support, c(8L, 9L))
  expect_identical(partition$residual_only_bucket_support, 7L)
  expect_identical(partition$summary$minimum_case_combo_count, 20L)
  expect_identical(partition$summary$residual_case_count, 11L)
  expect_equal(
    partition$summary$minimum_case_combo_share_of_non_near_threshold_misses,
    20 / 31,
    tolerance = 1e-12
  )
  expect_equal(
    partition$summary$residual_case_share_of_non_near_threshold_misses,
    11 / 31,
    tolerance = 1e-12
  )
  expect_true(isTRUE(
    partition$summary$combined_partition_reconstructs_non_near_threshold_misses
  ))
  expect_true(isTRUE(
    partition$summary$residual_frontier_starts_after_threshold
  ))
})

test_that("D12 partition-exclusive bucket helper stays executable in the seasonal-WCB package slice", {
  exclusive <- lwdid:::.story_local_wcb_partition_exclusive_bucket_representatives()

  expect_identical(exclusive$case_id, "TC-9.4.19")
  expect_identical(
    exclusive$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(exclusive$threshold_only_bucket_support, c(8L, 9L))
  expect_identical(exclusive$residual_only_bucket_support, 7L)
  expect_identical(exclusive$threshold_only_replay_seeds, c(72L, 74L))
  expect_identical(exclusive$residual_only_replay_seeds, c(48L, 53L, 79L))
  expect_identical(exclusive$summary$threshold_only_case_count, 2L)
  expect_identical(exclusive$summary$residual_only_case_count, 3L)
  expect_equal(
    exclusive$summary$threshold_only_share_of_partition_exclusive_cases,
    2 / 5,
    tolerance = 1e-12
  )
  expect_equal(
    exclusive$summary$residual_only_share_of_partition_exclusive_cases,
    3 / 5,
    tolerance = 1e-12
  )
  expect_identical(exclusive$summary$closest_threshold_only_replay_seed, 72L)
  expect_identical(exclusive$summary$highest_threshold_only_t_stat_seed, 74L)
  expect_identical(exclusive$summary$easiest_residual_only_replay_seed, 48L)
  expect_identical(exclusive$summary$hardest_residual_only_replay_seed, 79L)
  expect_true(isTRUE(exclusive$summary$all_partition_exclusive_targeted_replays_remain_misses))
  expect_true(isTRUE(exclusive$summary$bucket8_is_closest_threshold_only_case))
  expect_true(isTRUE(exclusive$summary$bucket9_has_highest_threshold_only_t_stat))
  expect_true(isTRUE(exclusive$summary$bucket7_reconstructs_residual_only_family))
})

test_that("D12 partition-exclusive signal contrast stays executable in the seasonal-WCB package slice", {
  contrast <- lwdid:::.story_local_wcb_partition_exclusive_bucket_signal_contrast()

  expect_identical(contrast$case_id, "TC-9.4.19")
  expect_identical(
    contrast$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(contrast$threshold_only_replay_seeds, c(72L, 74L))
  expect_identical(contrast$residual_only_replay_seeds, c(48L, 53L, 79L))
  expect_identical(
    contrast$threshold_only_pair$closest_threshold_only$replay_seed,
    72L
  )
  expect_identical(
    contrast$threshold_only_pair$highest_t_threshold_only$replay_seed,
    74L
  )
  expect_identical(
    contrast$residual_only_pair$easiest_residual_only$replay_seed,
    48L
  )
  expect_identical(
    contrast$residual_only_pair$hardest_residual_only$replay_seed,
    79L
  )
  expect_equal(
    contrast$summary$bucket9_t_stat_gain_over_bucket8,
    4.590259270017282,
    tolerance = 1e-15
  )
  expect_equal(
    contrast$summary$bucket8_t_stat_gain_over_easiest_residual,
    1.110600598370500,
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
})

test_that("D12 residual-only monotone decay stays executable in the seasonal-WCB package slice", {
  decay <- lwdid:::.story_local_wcb_partition_exclusive_residual_monotone_decay()

  expect_identical(decay$case_id, "TC-9.4.19")
  expect_identical(
    decay$story_live_blockers,
    c(
      "story-local-hardening-power-gap-tc-9-4-19",
      "d13-skip-backed-optional-backend-watchpoint"
    )
  )
  expect_identical(decay$threshold_only_replay_seeds, c(72L, 74L))
  expect_identical(decay$residual_only_replay_seeds, c(48L, 53L, 79L))
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
})

test_that("D12 residual-only targeted replay rank stability stays executable in the seasonal-WCB package slice", {
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
  expect_identical(
    rank_stability$replay_seed_offsets,
    c(6L, 11L, 37L)
  )
  expect_identical(
    rank_stability$replay_attempt_ids,
    c(7L, 12L, 38L)
  )
  expect_identical(
    rank_stability$replay_case_ranks,
    c(22L, 30L, 31L)
  )
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
  expect_equal(
    rank_stability$triplet_tuples,
    list(
      list(
        replay_seed = 48L,
        replay_attempt_id = 7L,
        replay_case_rank = 22L,
        replay_seed_offset = 6L,
        replication_id = 7L,
        treated_cluster_count = 7L,
        reject = FALSE,
        pvalue = 0.48828125,
        abs_gap_to_threshold = 0.43828125,
        abs_att = 0.482970992563102,
        original_se = 0.632549650013912,
        t_stat = 0.763530566418746
      ),
      list(
        replay_seed = 53L,
        replay_attempt_id = 12L,
        replay_case_rank = 30L,
        replay_seed_offset = 11L,
        replication_id = 12L,
        treated_cluster_count = 7L,
        reject = FALSE,
        pvalue = 0.8935546875,
        abs_gap_to_threshold = 0.8435546875,
        abs_att = 0.144728810559076,
        original_se = 0.909524337506003,
        t_stat = 0.159125824995443
      ),
      list(
        replay_seed = 79L,
        replay_attempt_id = 38L,
        replay_case_rank = 31L,
        replay_seed_offset = 37L,
        replication_id = 38L,
        treated_cluster_count = 7L,
        reject = FALSE,
        pvalue = 0.94921875,
        abs_gap_to_threshold = 0.89921875,
        abs_att = 0.0779617953550249,
        original_se = 1.12198708000843,
        t_stat = 0.0694854662269725
      )
    ),
    tolerance = 1e-12
  )
  expect_true(isTRUE(rank_stability$summary$replay_attempt_ids_strictly_ascend))
  expect_true(isTRUE(rank_stability$summary$replay_seed_offsets_strictly_ascend))
  expect_true(isTRUE(rank_stability$summary$case_ranks_strictly_ascend))
  expect_true(isTRUE(rank_stability$summary$abs_att_strictly_descends))
  expect_true(isTRUE(rank_stability$summary$original_se_strictly_ascends))
  expect_true(isTRUE(rank_stability$summary$t_stat_strictly_descends))
  expect_true(isTRUE(rank_stability$summary$pvalue_strictly_ascends))
  expect_true(isTRUE(rank_stability$summary$pvalue_gap_strictly_ascends))
  expect_true(isTRUE(rank_stability$summary$triplet_matches_residual_monotone_helper))
})

test_that("D1 rademacher bootstrap weights stay on {-1, +1}", {
  set.seed(42)
  weights <- lwdid:::.generate_weights(n_clusters = 256L, weight_type = "rademacher")

  expect_length(weights, 256L)
  expect_equal(sort(unique(weights)), c(-1, 1))
})

test_that("D2 Mammen bootstrap weights preserve the first three moments", {
  set.seed(42)
  weights <- lwdid:::.generate_weights(n_clusters = 50000L, weight_type = "mammen")

  expect_equal(mean(weights), 0, tolerance = 0.02)
  expect_equal(mean(weights^2), 1, tolerance = 0.02)
  expect_equal(mean(weights^3), 1, tolerance = 0.05)
})

test_that("D3 Webb bootstrap weights stay on the six-point support", {
  webb_support <- round(c(-sqrt(3 / 2), -1, -sqrt(1 / 2), sqrt(1 / 2), 1, sqrt(3 / 2)), 12)

  set.seed(42)
  weights <- lwdid:::.generate_weights(n_clusters = 50000L, weight_type = "webb")

  expect_equal(sort(unique(round(weights, 12))), webb_support)
})

test_that("D4 G=8 Rademacher full enumeration expands to 256 sign patterns", {
  weight_info <- lwdid:::.resolve_wild_bootstrap_weights(
    n_clusters = 8L,
    requested_n_bootstrap = 999L,
    weight_type = "rademacher",
    full_enumeration = NULL
  )

  expect_true(isTRUE(weight_info$full_enumeration))
  expect_identical(dim(weight_info$weights), c(256L, 8L))
  expect_identical(weight_info$actual_n_bootstrap, 256L)
  expect_equal(sort(unique(as.vector(weight_info$weights))), c(-1, 1))
})

test_that("D4 public native WCB auto-promotes the G=8 panel to 256 exact bootstrap draws", {
  expect_true(
    exists("make_e906_public_full_enumeration_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_full_enumeration_panel"
  )

  if (!exists("make_e906_public_full_enumeration_panel", mode = "function")) {
    return(invisible(NULL))
  }

  public_panel <- make_e906_public_full_enumeration_panel()
  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = public_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 1L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_true(isTRUE(wcb_result$full_enumeration))
  expect_identical(wcb_result$requested_n_bootstrap, 99L)
  expect_identical(wcb_result$actual_n_bootstrap, 256L)
  expect_length(wcb_result$att_bootstrap, 256L)
  expect_length(wcb_result$t_stats_bootstrap, 256L)
})

test_that("D5 full enumeration stays deterministic across seed changes", {
  set.seed(1)
  first_draw <- lwdid:::.resolve_wild_bootstrap_weights(
    n_clusters = 8L,
    requested_n_bootstrap = 999L,
    weight_type = "rademacher",
    full_enumeration = TRUE
  )

  set.seed(999)
  second_draw <- lwdid:::.resolve_wild_bootstrap_weights(
    n_clusters = 8L,
    requested_n_bootstrap = 999L,
    weight_type = "rademacher",
    full_enumeration = TRUE
  )

  expect_identical(first_draw$weights, second_draw$weights)
})

test_that("D5 public native WCB full enumeration is seed-invariant on the same G=8 panel", {
  expect_true(
    exists("make_e906_public_full_enumeration_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_full_enumeration_panel"
  )

  if (!exists("make_e906_public_full_enumeration_panel", mode = "function")) {
    return(invisible(NULL))
  }

  public_panel <- make_e906_public_full_enumeration_panel()
  first_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = public_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 1L,
      use_fwildclusterboot = FALSE
    )
  )
  second_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = public_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 999L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_true(isTRUE(first_result$full_enumeration))
  expect_true(isTRUE(second_result$full_enumeration))
  expect_identical(first_result$actual_n_bootstrap, 256L)
  expect_identical(second_result$actual_n_bootstrap, 256L)
  expect_identical(first_result$att_bootstrap, second_result$att_bootstrap)
  expect_identical(first_result$t_stats_bootstrap, second_result$t_stats_bootstrap)
  expect_identical(first_result$pvalue, second_result$pvalue)
  expect_identical(first_result$ci_lower, second_result$ci_lower)
  expect_identical(first_result$ci_upper, second_result$ci_upper)
})

test_that("D6 impose_null modes produce different bootstrap p-values on the same clustered panel", {
  core_panel <- make_e906_core_wcb_panel()
  impose_null_true <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = core_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 1L,
      full_enumeration = TRUE,
      use_fwildclusterboot = FALSE
    )
  )
  impose_null_false <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = core_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = FALSE,
      alpha = 0.05,
      seed = 1L,
      full_enumeration = TRUE,
      use_fwildclusterboot = FALSE
    )
  )

  expect_identical(impose_null_true$ci_method, "percentile_t")
  expect_identical(impose_null_false$ci_method, "percentile")
  expect_true(is.finite(impose_null_true$pvalue))
  expect_true(is.finite(impose_null_false$pvalue))
  expect_false(isTRUE(all.equal(impose_null_true$pvalue, impose_null_false$pvalue)))
})

test_that("D7 native WCB p-value uses the >= extremeness rule and stays inside [0, 1]", {
  core_panel <- make_e906_core_wcb_panel()
  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = core_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 1L,
      full_enumeration = TRUE,
      use_fwildclusterboot = FALSE
    )
  )

  valid_t <- wcb_result$t_stats_bootstrap[is.finite(wcb_result$t_stats_bootstrap)]
  pvalue_ge <- mean(abs(valid_t) >= abs(wcb_result$t_stat_original))
  pvalue_gt <- mean(abs(valid_t) > abs(wcb_result$t_stat_original))

  expect_gte(sum(abs(valid_t) == abs(wcb_result$t_stat_original)), 1L)
  expect_true(is.finite(wcb_result$pvalue))
  expect_gte(wcb_result$pvalue, 0)
  expect_lte(wcb_result$pvalue, 1)
  expect_equal(wcb_result$pvalue, pvalue_ge, tolerance = 1e-12)
  expect_gt(pvalue_ge, pvalue_gt)
})

test_that("D8 constant-outcome native WCB returns degenerate NaN inference surfaces", {
  constant_panel <- make_e906_constant_outcome_panel()
  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = constant_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 49L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(wcb_result, "lwdid_wcb_result")
  expect_lt(abs(wcb_result$original_se), 1e-10)
  expect_true(is.nan(wcb_result$pvalue))
  expect_true(is.nan(wcb_result$rejection_rate))
  expect_true(is.nan(wcb_result$ci_lower))
  expect_true(is.nan(wcb_result$ci_upper))
})

test_that("D9 native WCB original_se keeps the clustered small-sample correction", {
  core_panel <- make_e906_core_wcb_panel()
  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = core_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 1L,
      full_enumeration = TRUE,
      use_fwildclusterboot = FALSE
    )
  )

  X <- cbind("(Intercept)" = 1, D = core_panel$D)
  y <- core_panel$Y
  xtx_inv <- solve(crossprod(X))
  beta_hat <- drop(xtx_inv %*% crossprod(X, y))
  residuals <- y - drop(X %*% beta_hat)
  cluster_rows <- split(seq_len(nrow(core_panel)), core_panel$cluster)
  meat <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  for (rows in cluster_rows) {
    score_g <- drop(crossprod(X[rows, , drop = FALSE], residuals[rows]))
    meat <- meat + tcrossprod(score_g)
  }
  G <- length(cluster_rows)
  N <- nrow(core_panel)
  K <- ncol(X)
  correction <- (G / (G - 1)) * ((N - 1) / (N - K))
  expected_var <- correction * xtx_inv %*% meat %*% xtx_inv
  expected_se <- sqrt(max(as.numeric(expected_var[2, 2]), 0))

  expect_equal(wcb_result$original_se, expected_se, tolerance = 1e-12)
})

test_that("D14 restricted_model changes the imposed-null controls path when controls are present", {
  restricted_model_panel <- make_e906_restricted_model_panel()

  with_controls_result <- lwdid::wild_cluster_bootstrap(
    data = restricted_model_panel,
    y_transformed = "y",
    d = "D",
    cluster_var = "cluster",
    controls = c("control_a", "control_b"),
    n_bootstrap = 199L,
    weight_type = "rademacher",
    impose_null = TRUE,
    alpha = 0.05,
    seed = 1L,
    use_fwildclusterboot = FALSE,
    restricted_model = "with_controls"
  )
  intercept_only_result <- lwdid::wild_cluster_bootstrap(
    data = restricted_model_panel,
    y_transformed = "y",
    d = "D",
    cluster_var = "cluster",
    controls = c("control_a", "control_b"),
    n_bootstrap = 199L,
    weight_type = "rademacher",
    impose_null = TRUE,
    alpha = 0.05,
    seed = 1L,
    use_fwildclusterboot = FALSE,
    restricted_model = "intercept_only"
  )

  expect_s3_class(with_controls_result, "lwdid_wcb_result")
  expect_s3_class(intercept_only_result, "lwdid_wcb_result")
  expect_identical(with_controls_result$actual_n_bootstrap, 256L)
  expect_identical(intercept_only_result$actual_n_bootstrap, 256L)
  expect_identical(with_controls_result$restricted_model, "with_controls")
  expect_identical(intercept_only_result$restricted_model, "intercept_only")
  expect_equal(with_controls_result$att, intercept_only_result$att, tolerance = 1e-10)
  expect_equal(
    with_controls_result$original_se,
    intercept_only_result$original_se,
    tolerance = 1e-10
  )
  expect_identical(
    names(with_controls_result$restricted_coefficients),
    c("(Intercept)", "control_a", "control_b")
  )
  expect_identical(
    names(intercept_only_result$restricted_coefficients),
    "(Intercept)"
  )
  expect_false(isTRUE(all.equal(
    with_controls_result$fitted_base,
    intercept_only_result$fitted_base
  )))
  expect_false(isTRUE(all.equal(
    with_controls_result$resid_base,
    intercept_only_result$resid_base
  )))
  expect_true(is.finite(with_controls_result$pvalue))
  expect_true(is.finite(intercept_only_result$pvalue))
  expect_true(is.finite(with_controls_result$ci_lower))
  expect_true(is.finite(intercept_only_result$ci_lower))
  expect_false(isTRUE(all.equal(
    with_controls_result$pvalue,
    intercept_only_result$pvalue
  )))
  expect_false(isTRUE(all.equal(
    with_controls_result$ci_lower,
    intercept_only_result$ci_lower
  )))
  expect_false(isTRUE(all.equal(
    with_controls_result$ci_upper,
    intercept_only_result$ci_upper
  )))
})

test_that("D15 collinear WCB controls warn without aborting the native bootstrap path", {
  collinear_panel <- make_e906_collinear_panel()
  collinear_panel$cluster <- rep(c(1L, 1L, 2L, 2L), each = 10L)
  wcb_result <- NULL

  expect_warning(
    wcb_result <- lwdid::wild_cluster_bootstrap(
      data = collinear_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = c("control_a", "control_b"),
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 1L,
      use_fwildclusterboot = FALSE
    ),
    "condition number is large"
  )

  expect_s3_class(wcb_result, "lwdid_wcb_result")
  expect_true(is.finite(wcb_result$att))
  expect_true(is.finite(wcb_result$pvalue))
})

test_that("Task E9-06.7 singular native WCB controls keep a finite fallback path", {
  singular_panel <- make_e906_collinear_panel()
  singular_panel$control_c <- singular_panel$control_a
  singular_panel$cluster <- rep(c(1L, 1L, 2L, 2L), each = 10L)
  singular_result <- NULL

  expect_warning(
    singular_result <- lwdid::wild_cluster_bootstrap(
      data = singular_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = c("control_a", "control_b", "control_c"),
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 1L,
      use_fwildclusterboot = FALSE,
      restricted_model = "with_controls"
    ),
    "condition number is large"
  )

  expect_s3_class(singular_result, "lwdid_wcb_result")
  expect_identical(singular_result$actual_n_bootstrap, 4L)
  expect_identical(singular_result$restricted_model, "with_controls")
  expect_true(is.finite(singular_result$att))
  expect_true(is.finite(singular_result$original_se))
  expect_true(is.finite(singular_result$pvalue))
  expect_true(is.finite(singular_result$ci_lower))
  expect_true(is.finite(singular_result$ci_upper))
  expect_true(singular_result$ci_lower < singular_result$ci_upper)
  expect_true(all(is.finite(singular_result$fitted_base)))
  expect_true(all(is.finite(singular_result$resid_base)))
})

test_that("D18 native WCB keeps rejection_rate equal to pvalue after collinearity fallback", {
  collinear_panel <- make_e906_collinear_panel()
  collinear_panel$cluster <- rep(c(1L, 1L, 2L, 2L), each = 10L)

  wcb_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = collinear_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = c("control_a", "control_b"),
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 1L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(wcb_result, "lwdid_wcb_result")
  expect_equal(wcb_result$rejection_rate, wcb_result$pvalue, tolerance = 1e-12)
})

test_that("Task E9-06.0 fixtures feed the public wild_cluster_bootstrap surface", {
  public_result <- run_seasonal_wcb_public_smoke(seed = 99L)

  expect_s3_class(public_result, "lwdid_wcb_result")
  expect_identical(public_result$requested_n_bootstrap, 199L)
  expect_identical(public_result$actual_n_bootstrap, 1024L)
  expect_true(isTRUE(public_result$full_enumeration))
  expect_true(is.finite(public_result$att))
  expect_true(is.finite(public_result$pvalue))
  expect_lt(public_result$ci_lower, public_result$ci_upper)
})

test_that("E1 public lwdid demeanq returns a finite ATT on a quarterly seasonal panel", {
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  panel_data <- make_e906_public_lwdid_panel()
  demeanq_result <- suppressWarnings(
    lwdid::lwdid(
      data = panel_data,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "demeanq",
      season_var = "quarter",
      Q = 4L
    )
  )

  expect_s3_class(demeanq_result, "lwdid_result")
  expect_identical(demeanq_result$rolling, "demeanq")
  expect_identical(demeanq_result$exclude_pre_periods, 0L)
  expect_true(is.finite(demeanq_result$att))
  expect_equal(demeanq_result$att, 1.25, tolerance = 1e-4)
})

test_that("E2 public lwdid detrendq returns a finite ATT on the same quarterly panel", {
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  panel_data <- make_e906_public_lwdid_panel()
  detrendq_result <- suppressWarnings(
    lwdid::lwdid(
      data = panel_data,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "detrendq",
      season_var = "quarter",
      Q = 4L
    )
  )

  expect_s3_class(detrendq_result, "lwdid_result")
  expect_identical(detrendq_result$rolling, "detrendq")
  expect_identical(detrendq_result$exclude_pre_periods, 0L)
  expect_true(is.finite(detrendq_result$att))
  expect_equal(detrendq_result$att, 1.250375, tolerance = 1e-4)
  expect_true(is.finite(detrendq_result$pvalue))
})

test_that("E3 public lwdid WCB returns finite p-value and CI on the clustered panel", {
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  panel_data <- make_e906_public_lwdid_panel()
  expect_true(
    "cluster" %in% names(panel_data),
    info = "public WCB integration panel must expose a cluster column"
  )

  if (!("cluster" %in% names(panel_data))) {
    return(invisible(NULL))
  }

  wcb_result <- suppressWarnings(
    lwdid::lwdid(
      data = panel_data,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      vce = "bootstrap",
      cluster_var = "cluster",
      wcb_reps = 19L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(wcb_result, "lwdid_result")
  expect_identical(wcb_result$vce_type, "wild_cluster_bootstrap")
  expect_true(is.finite(wcb_result$att))
  expect_true(is.finite(wcb_result$pvalue))
  expect_true(is.finite(wcb_result$ci_lower))
  expect_true(is.finite(wcb_result$ci_upper))
  expect_lt(wcb_result$ci_lower, wcb_result$att)
  expect_gt(wcb_result$ci_upper, wcb_result$att)
  expect_s3_class(wcb_result$wcb_details, "lwdid_wcb_result")
  expect_identical(wcb_result$wcb_details$requested_n_bootstrap, 19L)
})

test_that("E4 public lwdid demeanq plus WCB stays finite on the clustered seasonal panel", {
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  panel_data <- make_e906_public_lwdid_panel()
  expect_true(
    "cluster" %in% names(panel_data),
    info = "public clustered seasonal panel must expose a cluster column"
  )

  if (!("cluster" %in% names(panel_data))) {
    return(invisible(NULL))
  }

  demeanq_wcb_result <- suppressWarnings(
    lwdid::lwdid(
      data = panel_data,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "demeanq",
      season_var = "quarter",
      Q = 4L,
      vce = "bootstrap",
      cluster_var = "cluster",
      wcb_reps = 19L,
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(demeanq_wcb_result, "lwdid_result")
  expect_identical(demeanq_wcb_result$rolling, "demeanq")
  expect_identical(demeanq_wcb_result$vce_type, "wild_cluster_bootstrap")
  expect_true(is.finite(demeanq_wcb_result$att))
  expect_equal(demeanq_wcb_result$att, 1.25, tolerance = 1e-4)
  expect_true(is.finite(demeanq_wcb_result$pvalue))
  expect_true(is.finite(demeanq_wcb_result$ci_lower))
  expect_true(is.finite(demeanq_wcb_result$ci_upper))
  expect_lt(demeanq_wcb_result$ci_lower, demeanq_wcb_result$ci_upper)
  expect_s3_class(demeanq_wcb_result$wcb_details, "lwdid_wcb_result")
})

test_that("E5 clustered public lwdid auto-triggers WCB when G is below 20", {
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  panel_data <- make_e906_public_lwdid_panel()
  expect_true(
    "cluster" %in% names(panel_data),
    info = "public auto-WCB panel must expose a cluster column"
  )

  if (!("cluster" %in% names(panel_data))) {
    return(invisible(NULL))
  }

  expect_lt(length(unique(panel_data$cluster)), 20L)

  expect_message(
    auto_wcb_result <- suppressWarnings(
      lwdid::lwdid(
        data = panel_data,
        y = "y",
        ivar = "id",
        tvar = "time",
        d = "treat",
        post = "post",
        vce = "cluster",
        cluster_var = "cluster",
        auto_wcb = TRUE,
        wcb_reps = 19L,
        use_fwildclusterboot = FALSE
      )
    ),
    regexp = "Wild Cluster Bootstrap"
  )

  expect_s3_class(auto_wcb_result, "lwdid_result")
  expect_true(isTRUE(auto_wcb_result$wcb_auto_triggered))
  expect_identical(auto_wcb_result$vce_type, "wild_cluster_bootstrap")
  expect_true(is.finite(auto_wcb_result$pvalue))
  expect_true(is.finite(auto_wcb_result$ci_lower))
  expect_true(is.finite(auto_wcb_result$ci_upper))
  expect_s3_class(auto_wcb_result$wcb_details, "lwdid_wcb_result")
})

test_that("E6 public lwdid exclude_pre_periods removes the contaminated late pre-period seasonal fit", {
  expect_true(
    exists("make_e906_public_lwdid_anticipation_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_anticipation_panel"
  )

  if (!exists("make_e906_public_lwdid_anticipation_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contaminated_panel <- make_e906_public_lwdid_anticipation_panel()
  baseline_result <- suppressWarnings(
    lwdid::lwdid(
      data = contaminated_panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "demeanq",
      season_var = "quarter",
      Q = 4L
    )
  )
  excluded_result <- suppressWarnings(
    lwdid::lwdid(
      data = contaminated_panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "demeanq",
      season_var = "quarter",
      Q = 4L,
      exclude_pre_periods = 2L
    )
  )

  expect_s3_class(baseline_result, "lwdid_result")
  expect_s3_class(excluded_result, "lwdid_result")
  expect_identical(excluded_result$exclude_pre_periods, 2L)
  expect_true(is.finite(baseline_result$att))
  expect_true(is.finite(excluded_result$att))
  expect_gt(abs(baseline_result$att - excluded_result$att), 0.15)
  expect_equal(excluded_result$att, 1.25, tolerance = 1e-4)
  expect_true(is.finite(excluded_result$pvalue))
})

test_that("Task E9-06.13 quarter alias WCB keeps the clustered seasonal public surface stable", {
  expect_true(
    exists("make_e906_post_consume_closure_wave_public_contracts", mode = "function"),
    info = "missing Task E9-06.13 helper: make_e906_post_consume_closure_wave_public_contracts"
  )
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_post_consume_closure_wave_public_contracts", mode = "function") ||
    !exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_post_consume_closure_wave_public_contracts()
  panel_data <- make_e906_public_lwdid_panel()
  result <- suppressWarnings(
    lwdid::lwdid(
      data = panel_data,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "demeanq",
      quarter = "quarter",
      Q = 4L,
      vce = "bootstrap",
      cluster_var = "cluster",
      wcb_reps = 19L,
      use_fwildclusterboot = FALSE
    )
  )

  expected <- contracts$quarter_alias_demeanq_wcb
  expect_s3_class(result, "lwdid_result")
  expect_identical(result$rolling, "demeanq")
  expect_identical(result$vce_type, "wild_cluster_bootstrap")
  expect_equal(result$att, expected$att, tolerance = 1e-12)
  expect_equal(result$pvalue, expected$pvalue, tolerance = 1e-12)
  expect_equal(result$ci_lower, expected$ci_lower, tolerance = 1e-12)
  expect_equal(result$ci_upper, expected$ci_upper, tolerance = 1e-12)
  expect_identical(result$wcb_details$actual_n_bootstrap, expected$actual_n_bootstrap)
  expect_identical(result$wcb_details$requested_n_bootstrap, expected$requested_n_bootstrap)
})

test_that("Task E9-06.13 public detrendq plus WCB stays finite on the clustered seasonal panel", {
  expect_true(
    exists("make_e906_post_consume_closure_wave_public_contracts", mode = "function"),
    info = "missing Task E9-06.13 helper: make_e906_post_consume_closure_wave_public_contracts"
  )
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_post_consume_closure_wave_public_contracts", mode = "function") ||
    !exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_post_consume_closure_wave_public_contracts()
  panel_data <- make_e906_public_lwdid_panel()
  result <- suppressWarnings(
    lwdid::lwdid(
      data = panel_data,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "detrendq",
      season_var = "quarter",
      Q = 4L,
      vce = "bootstrap",
      cluster_var = "cluster",
      wcb_reps = 19L,
      use_fwildclusterboot = FALSE
    )
  )

  expected <- contracts$detrendq_wcb
  expect_s3_class(result, "lwdid_result")
  expect_identical(result$rolling, "detrendq")
  expect_identical(result$vce_type, "wild_cluster_bootstrap")
  expect_equal(result$att, expected$att, tolerance = 1e-12)
  expect_equal(result$pvalue, expected$pvalue, tolerance = 1e-12)
  expect_equal(result$ci_lower, expected$ci_lower, tolerance = 1e-12)
  expect_equal(result$ci_upper, expected$ci_upper, tolerance = 1e-12)
  expect_identical(result$wcb_details$actual_n_bootstrap, expected$actual_n_bootstrap)
  expect_identical(result$wcb_details$requested_n_bootstrap, expected$requested_n_bootstrap)
})

test_that("Task E9-06.13 public detrendq auto-WCB honors exclude_pre_periods on the anticipation panel", {
  expect_true(
    exists("make_e906_post_consume_closure_wave_public_contracts", mode = "function"),
    info = "missing Task E9-06.13 helper: make_e906_post_consume_closure_wave_public_contracts"
  )
  expect_true(
    exists("make_e906_public_lwdid_anticipation_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_anticipation_panel"
  )

  if (!exists("make_e906_post_consume_closure_wave_public_contracts", mode = "function") ||
    !exists("make_e906_public_lwdid_anticipation_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_post_consume_closure_wave_public_contracts()
  contaminated_panel <- make_e906_public_lwdid_anticipation_panel()
  baseline_result <- suppressWarnings(
    lwdid::lwdid(
      data = contaminated_panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "detrendq",
      season_var = "quarter",
      Q = 4L,
      vce = "bootstrap",
      cluster_var = "cluster",
      wcb_reps = 19L,
      use_fwildclusterboot = FALSE
    )
  )
  excluded_result <- suppressWarnings(
    lwdid::lwdid(
      data = contaminated_panel,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "treat",
      post = "post",
      rolling = "detrendq",
      season_var = "quarter",
      Q = 4L,
      exclude_pre_periods = 2L,
      vce = "cluster",
      cluster_var = "cluster",
      auto_wcb = TRUE,
      wcb_reps = 19L,
      use_fwildclusterboot = FALSE
    )
  )

  expected_baseline <- contracts$detrendq_bootstrap_anticipation
  expected_excluded <- contracts$detrendq_excluded_auto_wcb
  expect_s3_class(baseline_result, "lwdid_result")
  expect_s3_class(excluded_result, "lwdid_result")
  expect_identical(excluded_result$rolling, "detrendq")
  expect_identical(excluded_result$exclude_pre_periods, 2L)
  expect_true(isTRUE(excluded_result$wcb_auto_triggered))
  expect_equal(baseline_result$att, expected_baseline$att, tolerance = 1e-12)
  expect_equal(baseline_result$pvalue, expected_baseline$pvalue, tolerance = 1e-12)
  expect_equal(baseline_result$ci_lower, expected_baseline$ci_lower, tolerance = 1e-12)
  expect_equal(baseline_result$ci_upper, expected_baseline$ci_upper, tolerance = 1e-12)
  expect_equal(excluded_result$att, expected_excluded$att, tolerance = 1e-12)
  expect_equal(excluded_result$pvalue, expected_excluded$pvalue, tolerance = 1e-12)
  expect_equal(excluded_result$ci_lower, expected_excluded$ci_lower, tolerance = 1e-12)
  expect_equal(excluded_result$ci_upper, expected_excluded$ci_upper, tolerance = 1e-12)
  expect_identical(excluded_result$wcb_details$actual_n_bootstrap, expected_excluded$actual_n_bootstrap)
  expect_gt(excluded_result$att - baseline_result$att, 0.7)
})

test_that("Task E9-06.13.5 B5 monthly demeanq keeps the exact Q=12 package surface", {
  expect_true(
    exists("make_e906_forced_bundle_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.5 helper: make_e906_forced_bundle_package_contracts"
  )
  expect_true(
    exists("make_e906_monthly_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_monthly_panel"
  )

  if (!exists("make_e906_forced_bundle_package_contracts", mode = "function") ||
    !exists("make_e906_monthly_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_forced_bundle_package_contracts()
  monthly_panel <- make_e906_monthly_panel()
  monthly_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(monthly_panel),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 14L,
    rolling = "demeanq",
    post = "post",
    season_var = "month",
    Q = 12L
  )

  expected <- contracts$monthly_demeanq
  expect_equal(monthly_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(monthly_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(monthly_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(monthly_result$firstpost, expected$firstpost)
  expect_identical(monthly_result$n_pre, expected$n_pre)
})

test_that("Task E9-06.13.5 B3 sparse seasonal units keep the exact package-surface validation error", {
  expect_true(
    exists("make_e906_forced_bundle_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.5 helper: make_e906_forced_bundle_package_contracts"
  )
  expect_true(
    exists("make_e906_sparse_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_sparse_panel"
  )

  if (!exists("make_e906_forced_bundle_package_contracts", mode = "function") ||
    !exists("make_e906_sparse_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_forced_bundle_package_contracts()
  sparse_panel <- make_e906_sparse_panel()

  expect_error(
    lwdid:::transform_common(
      dt = data.table::as.data.table(sparse_panel),
      y = "y",
      ivar = "id",
      tvar = "tindex",
      g = 3L,
      rolling = "demeanq",
      post = "post",
      season_var = "quarter",
      Q = 4L
    ),
    contracts$sparse_demeanq$error,
    fixed = TRUE
  )
})

test_that("Task E9-06.13.5 C11 detrendq and C12 quarter alias keep the exact excluded package surface", {
  expect_true(
    exists("make_e906_forced_bundle_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.5 helper: make_e906_forced_bundle_package_contracts"
  )
  expect_true(
    exists("make_e906_panel_data", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_panel_data"
  )

  if (!exists("make_e906_forced_bundle_package_contracts", mode = "function") ||
    !exists("make_e906_panel_data", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_forced_bundle_package_contracts()
  panel_data <- make_e906_panel_data()
  alias_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    rolling = "detrendq",
    post = "post",
    quarter = "quarter",
    Q = 4L,
    exclude_pre_periods = 2L
  )

  expected <- contracts$detrendq_alias_excluded
  expect_equal(alias_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(alias_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(alias_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(alias_result$firstpost, expected$firstpost)
  expect_identical(alias_result$n_pre, expected$n_pre)
})

test_that("Task E9-06.13.5 F3 biannual demeanq keeps exact post means and firstpost on the minimal-Q package surface", {
  expect_true(
    exists("make_e906_forced_bundle_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.5 helper: make_e906_forced_bundle_package_contracts"
  )
  expect_true(
    exists("make_e906_biannual_identity_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_biannual_identity_panel"
  )

  if (!exists("make_e906_forced_bundle_package_contracts", mode = "function") ||
    !exists("make_e906_biannual_identity_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_forced_bundle_package_contracts()
  biannual_panel <- make_e906_biannual_identity_panel()
  biannual_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(biannual_panel),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 7L,
    rolling = "demeanq",
    post = "post",
    season_var = "halfyear",
    Q = 2L
  )

  expected <- contracts$biannual_demeanq
  expect_equal(biannual_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(biannual_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(biannual_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(biannual_result$firstpost, expected$firstpost)
  expect_identical(biannual_result$n_pre, expected$n_pre)
})

test_that("Task E9-06.13.6 Task E9-06.7 demeanq keeps the exact post-summary package surface", {
  expect_true(
    exists("make_e906_post_summary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.6 helper: make_e906_post_summary_package_contracts"
  )
  expect_true(
    exists("make_e906_panel_data", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_panel_data"
  )

  if (!exists("make_e906_post_summary_package_contracts", mode = "function") ||
    !exists("make_e906_panel_data", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_post_summary_package_contracts()
  panel_data <- make_e906_panel_data()
  demeanq_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    rolling = "demeanq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )

  expected <- contracts$demeanq_base
  expect_match(expected$source_anchor, "Task E9-06.7", fixed = TRUE)
  expect_equal(demeanq_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(demeanq_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(demeanq_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(demeanq_result$firstpost, expected$firstpost)
  expect_identical(demeanq_result$n_pre, expected$n_pre)
})

test_that("Task E9-06.13.6 B6 demeanq exclude_pre_periods keeps the exact shifted package surface", {
  expect_true(
    exists("make_e906_post_summary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.6 helper: make_e906_post_summary_package_contracts"
  )
  expect_true(
    exists("make_e906_panel_data", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_panel_data"
  )

  if (!exists("make_e906_post_summary_package_contracts", mode = "function") ||
    !exists("make_e906_panel_data", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_post_summary_package_contracts()
  panel_data <- make_e906_panel_data()
  demeanq_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    rolling = "demeanq",
    post = "post",
    season_var = "quarter",
    Q = 4L,
    exclude_pre_periods = 2L
  )

  expected <- contracts$demeanq_excluded
  expect_match(expected$source_anchor, "B6 exclude_pre_periods", fixed = TRUE)
  expect_equal(demeanq_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(demeanq_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(demeanq_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(demeanq_result$firstpost, expected$firstpost)
  expect_identical(demeanq_result$n_pre, expected$n_pre)
})

test_that("Task E9-06.13.6 C8/C9 detrendq keeps the exact post-summary package surface", {
  expect_true(
    exists("make_e906_post_summary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.6 helper: make_e906_post_summary_package_contracts"
  )
  expect_true(
    exists("make_e906_panel_data", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_panel_data"
  )

  if (!exists("make_e906_post_summary_package_contracts", mode = "function") ||
    !exists("make_e906_panel_data", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_post_summary_package_contracts()
  panel_data <- make_e906_panel_data()
  detrendq_result <- lwdid:::transform_common(
    dt = data.table::as.data.table(panel_data),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 9L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )

  expected <- contracts$detrendq_base
  expect_match(expected$source_anchor, "C8/C9", fixed = TRUE)
  expect_equal(detrendq_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(detrendq_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(detrendq_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(detrendq_result$firstpost, expected$firstpost)
  expect_identical(detrendq_result$n_pre, expected$n_pre)
})

test_that("Task E9-06.13.7 A4 auto-detect keeps the exact seasonal frequency package surface", {
  expect_true(
    exists("make_e906_weekly_frequency_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.7 helper: make_e906_weekly_frequency_package_contracts"
  )
  expect_true(
    exists("make_e906_frequency_detection_cases", mode = "function"),
    info = "missing Task E9-06.1 detection-case helper"
  )

  if (!exists("make_e906_weekly_frequency_package_contracts", mode = "function") ||
    !exists("make_e906_frequency_detection_cases", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_weekly_frequency_package_contracts()
  detection_cases <- make_e906_frequency_detection_cases()

  for (case_name in names(contracts$frequency_detection)) {
    detection <- lwdid:::.auto_detect_frequency(
      data = detection_cases[[case_name]]$data,
      tvar = "year",
      ivar = "id"
    )
    expected <- contracts$frequency_detection[[case_name]]

    expect_identical(detection$frequency, expected$frequency, info = case_name)
    expect_identical(detection$Q, expected$Q, info = case_name)
    expect_identical(detection$confidence, expected$confidence, info = case_name)
  }
})

test_that("Task E9-06.13.7 C7 weekly detrend keeps the exact package surface", {
  expect_true(
    exists("make_e906_weekly_frequency_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.7 helper: make_e906_weekly_frequency_package_contracts"
  )
  expect_true(
    exists("make_e906_weekly_detrend_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_weekly_detrend_panel"
  )

  if (!exists("make_e906_weekly_frequency_package_contracts", mode = "function") ||
    !exists("make_e906_weekly_detrend_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_weekly_frequency_package_contracts()
  weekly_panel <- make_e906_weekly_detrend_panel()
  weekly_result <- lwdid:::transform_detrendq(
    dt = data.table::as.data.table(weekly_panel),
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 55L,
    season_var = "week",
    Q = 52L,
    post = "post"
  )

  expected <- contracts$weekly_detrend
  expect_match(expected$source_anchor, "C7", fixed = TRUE)
  expect_equal(weekly_result$seasonal_fit, expected$seasonal_fit, tolerance = 1e-10)
  expect_equal(weekly_result$y_trans, expected$y_trans, tolerance = 1e-10)
  expect_equal(weekly_result$ydot_postavg, expected$ydot_postavg, tolerance = 1e-10)
  expect_identical(weekly_result$firstpost, expected$firstpost)
  expect_identical(weekly_result$n_pre, expected$n_pre)
})

test_that("Task E9-06.13.7 C10 invalid seasonal values keep the exact dispatch validation error", {
  expect_true(
    exists("make_e906_weekly_frequency_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.7 helper: make_e906_weekly_frequency_package_contracts"
  )
  expect_true(
    exists("make_e906_panel_data", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_panel_data"
  )

  if (!exists("make_e906_weekly_frequency_package_contracts", mode = "function") ||
    !exists("make_e906_panel_data", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_weekly_frequency_package_contracts()
  panel_data <- make_e906_panel_data()
  panel_data$quarter[[1L]] <- 5L

  expect_error(
    lwdid:::transform_common(
      dt = data.table::as.data.table(panel_data),
      y = "y",
      ivar = "id",
      tvar = "tindex",
      g = 9L,
      rolling = "detrendq",
      post = "post",
      season_var = "quarter",
      Q = 4L
    ),
    contracts$invalid_detrendq_error$message,
    fixed = TRUE
  )
})

test_that("Task E9-06.13.8 A5 keeps the exact Q-boundary validation error", {
  expect_true(
    exists("make_e906_validation_boundary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.8 helper: make_e906_validation_boundary_package_contracts"
  )

  if (!exists("make_e906_validation_boundary_package_contracts", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_validation_boundary_package_contracts()
  q_lt_two_panel <- data.frame(
    id = rep(1L, 4L),
    tindex = 1:4,
    quarter = rep(1L, 4L),
    post = c(0L, 0L, 1L, 1L),
    y = c(10.0, 10.5, 11.0, 11.5)
  )

  q_error <- tryCatch(
    {
      lwdid:::.validate_seasonal_inputs(
        data = q_lt_two_panel,
        y = "y",
        season_var = "quarter",
        Q = 1L,
        ivar = "id",
        tindex = "tindex",
        post = "post",
        method = "demeanq",
        exclude_pre_periods = 0L,
        min_global_pre_periods = 1L
      )
      NULL
    },
    error = function(err) err
  )

  expected <- contracts$q_lt_two_error
  expect_match(expected$source_anchor, "A5", fixed = TRUE)
  expect_s3_class(q_error, expected$class, exact = FALSE)
  expect_identical(conditionMessage(q_error), expected$message)
})

test_that("Task E9-06.13.8 A2 keeps the exact seasonal coverage warning surface", {
  expect_true(
    exists("make_e906_validation_boundary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.8 helper: make_e906_validation_boundary_package_contracts"
  )
  expect_true(
    exists("make_e906_validation_case_bundle", mode = "function"),
    info = "missing Task E9-06.1 validation-case helper"
  )

  if (!exists("make_e906_validation_boundary_package_contracts", mode = "function") ||
    !exists("make_e906_validation_case_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_validation_boundary_package_contracts()
  validation_cases <- make_e906_validation_case_bundle()
  warning_messages <- character(0)
  a2_result <- withCallingHandlers(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a2_coverage_warning$call),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expected <- contracts$coverage_warning
  expect_match(expected$source_anchor, "A2", fixed = TRUE)
  expect_identical(warning_messages, expected$warnings)
  expect_identical(a2_result$warnings, expected$warnings)
  expect_true(isTRUE(a2_result$is_valid))
  expect_identical(a2_result$K, expected$K)
  expect_identical(a2_result$n_units_valid, expected$n_units_valid)
  expect_identical(a2_result$n_units_invalid, expected$n_units_invalid)
  expect_identical(a2_result$season_coverage$all_covered, expected$all_covered)
  expect_identical(
    unname(a2_result$season_coverage$units_with_gaps[[1L]]),
    expected$units_with_gaps[[1L]]
  )
  expect_identical(a2_result$season_coverage$total_uncovered, expected$total_uncovered)
})

test_that("Task E9-06.13.8 A6/A7 keep the exact pre-window validation boundaries", {
  expect_true(
    exists("make_e906_validation_boundary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.8 helper: make_e906_validation_boundary_package_contracts"
  )
  expect_true(
    exists("make_e906_validation_case_bundle", mode = "function"),
    info = "missing Task E9-06.1 validation-case helper"
  )

  if (!exists("make_e906_validation_boundary_package_contracts", mode = "function") ||
    !exists("make_e906_validation_case_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_validation_boundary_package_contracts()
  validation_cases <- make_e906_validation_case_bundle()

  a6_base <- do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a6_exclude_pre_base$call)
  a6_shifted_error <- tryCatch(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a6_exclude_pre_shifted$call),
    error = function(err) err
  )
  a7_error <- tryCatch(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a7_min_global_pre$call),
    error = function(err) err
  )

  expected_base <- contracts$exclude_pre_base
  expect_match(expected_base$source_anchor, "A6", fixed = TRUE)
  expect_true(isTRUE(a6_base$is_valid))
  expect_identical(a6_base$K, expected_base$K)
  expect_identical(a6_base$n_units_valid, expected_base$n_units_valid)
  expect_identical(a6_base$n_units_invalid, expected_base$n_units_invalid)
  expect_identical(a6_base$season_coverage$all_covered, expected_base$all_covered)
  expect_identical(a6_base$season_coverage$total_uncovered, expected_base$total_uncovered)

  expected_shifted <- contracts$exclude_pre_shifted_error
  expect_s3_class(a6_shifted_error, expected_shifted$class, exact = FALSE)
  expect_identical(conditionMessage(a6_shifted_error), expected_shifted$message)

  expected_a7 <- contracts$min_global_pre_error
  expect_match(expected_a7$source_anchor, "A7", fixed = TRUE)
  expect_s3_class(a7_error, expected_a7$class, exact = FALSE)
  expect_identical(conditionMessage(a7_error), expected_a7$message)
})

test_that("Task E9-06.13.8 A8 keeps the exact quarter-alias validation surface", {
  expect_true(
    exists("make_e906_validation_boundary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.8 helper: make_e906_validation_boundary_package_contracts"
  )
  expect_true(
    exists("make_e906_validation_case_bundle", mode = "function"),
    info = "missing Task E9-06.1 validation-case helper"
  )

  if (!exists("make_e906_validation_boundary_package_contracts", mode = "function") ||
    !exists("make_e906_validation_case_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_validation_boundary_package_contracts()
  validation_cases <- make_e906_validation_case_bundle()
  a8_result <- do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a8_quarter_alias$call)

  expected <- contracts$quarter_alias
  expect_match(expected$source_anchor, "A8", fixed = TRUE)
  expect_true(isTRUE(a8_result$is_valid))
  expect_identical(a8_result$K, expected$K)
  expect_identical(a8_result$n_units_valid, expected$n_units_valid)
  expect_identical(a8_result$n_units_invalid, expected$n_units_invalid)
  expect_identical(a8_result$season_coverage$all_covered, expected$all_covered)
  expect_identical(a8_result$season_coverage$total_uncovered, expected$total_uncovered)
  expect_identical(a8_result$warnings, expected$warnings)
})

test_that("Task E9-06.13.9 A1 keeps the exact invalid-range validation error", {
  expect_true(
    exists("make_e906_validation_error_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.9 helper: make_e906_validation_error_package_contracts"
  )
  expect_true(
    exists("make_e906_validation_case_bundle", mode = "function"),
    info = "missing Task E9-06.1 validation-case helper"
  )

  if (!exists("make_e906_validation_error_package_contracts", mode = "function") ||
    !exists("make_e906_validation_case_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_validation_error_package_contracts()
  validation_cases <- make_e906_validation_case_bundle()
  a1_error <- tryCatch(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a1_invalid_range$call),
    error = function(err) err
  )

  expected <- contracts$invalid_range_error
  expect_match(expected$source_anchor, "A1", fixed = TRUE)
  expect_s3_class(a1_error, expected$class, exact = FALSE)
  expect_identical(conditionMessage(a1_error), expected$message)
})

test_that("Task E9-06.13.9 A3 keeps the exact insufficient-pre validation error", {
  expect_true(
    exists("make_e906_validation_error_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.9 helper: make_e906_validation_error_package_contracts"
  )
  expect_true(
    exists("make_e906_validation_case_bundle", mode = "function"),
    info = "missing Task E9-06.1 validation-case helper"
  )

  if (!exists("make_e906_validation_error_package_contracts", mode = "function") ||
    !exists("make_e906_validation_case_bundle", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_validation_error_package_contracts()
  validation_cases <- make_e906_validation_case_bundle()
  a3_error <- tryCatch(
    do.call(lwdid:::.validate_seasonal_inputs, validation_cases$a3_insufficient_pre$call),
    error = function(err) err
  )

  expected <- contracts$insufficient_pre_error
  expect_match(expected$source_anchor, "A3", fixed = TRUE)
  expect_s3_class(a3_error, expected$class, exact = FALSE)
  expect_identical(conditionMessage(a3_error), expected$message)
})

test_that("Task E9-06.13.10 D16 keeps the exact test_inversion CI package surface", {
  expect_true(
    exists("make_e906_test_inversion_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.10 helper: make_e906_test_inversion_package_contracts"
  )
  expect_true(
    exists("make_e906_test_inversion_shared_fixture", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_test_inversion_shared_fixture"
  )

  if (!exists("make_e906_test_inversion_package_contracts", mode = "function") ||
    !exists("make_e906_test_inversion_shared_fixture", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_test_inversion_package_contracts()
  fixture <- make_e906_test_inversion_shared_fixture()
  ti_result <- suppressWarnings(
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

  expected <- contracts$ci_surface
  expect_match(expected$source_anchor, "D16", fixed = TRUE)
  expect_identical(ti_result$ci_method, expected$ci_method)
  expect_identical(ti_result$actual_n_bootstrap, expected$actual_n_bootstrap)
  expect_equal(ti_result$att, expected$att, tolerance = 1e-12)
  expect_equal(ti_result$ci_lower, expected$ci_lower, tolerance = 1e-12)
  expect_equal(ti_result$ci_upper, expected$ci_upper, tolerance = 1e-12)
  expect_equal(ti_result$pvalue, expected$pvalue, tolerance = 1e-12)
  expect_lte(ti_result$ci_lower, ti_result$att)
  expect_gte(ti_result$ci_upper, ti_result$att)
})

test_that("Task E9-06.13.10 D17 keeps the exact test_inversion width-band package surface", {
  expect_true(
    exists("make_e906_test_inversion_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.10 helper: make_e906_test_inversion_package_contracts"
  )
  expect_true(
    exists("make_e906_test_inversion_shared_fixture", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_test_inversion_shared_fixture"
  )

  if (!exists("make_e906_test_inversion_package_contracts", mode = "function") ||
    !exists("make_e906_test_inversion_shared_fixture", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_test_inversion_package_contracts()
  fixture <- make_e906_test_inversion_shared_fixture()
  test_inversion <- suppressWarnings(
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
  percentile_t <- suppressWarnings(
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

  width_test_inversion <- test_inversion$ci_upper - test_inversion$ci_lower
  width_percentile_t <- percentile_t$ci_upper - percentile_t$ci_lower
  width_ratio <- width_test_inversion / width_percentile_t
  relative_gap <- abs(width_test_inversion - width_percentile_t) / abs(width_percentile_t)

  expected <- contracts$width_surface
  expect_match(expected$source_anchor, "D17", fixed = TRUE)
  expect_identical(test_inversion$actual_n_bootstrap, expected$test_inversion_actual_n_bootstrap)
  expect_identical(percentile_t$actual_n_bootstrap, expected$percentile_t_actual_n_bootstrap)
  expect_equal(width_ratio, expected$width_ratio, tolerance = 1e-12)
  expect_equal(relative_gap, expected$relative_gap, tolerance = 1e-12)
  expect_gte(width_ratio, 0.5)
  expect_lte(width_ratio, 1.5)
  expect_lte(relative_gap, 0.5)
})

test_that("Task E9-06.13.10 D19 keeps the exact test_inversion monotonicity package surface", {
  expect_true(
    exists("make_e906_test_inversion_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.10 helper: make_e906_test_inversion_package_contracts"
  )
  expect_true(
    exists("make_e906_test_inversion_shared_fixture", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_test_inversion_shared_fixture"
  )

  if (!exists("make_e906_test_inversion_package_contracts", mode = "function") ||
    !exists("make_e906_test_inversion_shared_fixture", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_test_inversion_package_contracts()
  fixture <- make_e906_test_inversion_shared_fixture()
  ti_result <- suppressWarnings(
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
  pvalue_at_att <- lwdid:::.compute_bootstrap_pvalue_at_null(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    null_value = ti_result$att,
    att_original = ti_result$att,
    se_original = ti_result$original_se,
    n_bootstrap = 999L,
    weight_type = "rademacher",
    seed = 123L
  )
  pvalue_far_upper <- lwdid:::.compute_bootstrap_pvalue_at_null(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    null_value = ti_result$att + 5 * ti_result$original_se,
    att_original = ti_result$att,
    se_original = ti_result$original_se,
    n_bootstrap = 999L,
    weight_type = "rademacher",
    seed = 123L
  )
  pvalue_far_lower <- lwdid:::.compute_bootstrap_pvalue_at_null(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    null_value = ti_result$att - 5 * ti_result$original_se,
    att_original = ti_result$att,
    se_original = ti_result$original_se,
    n_bootstrap = 999L,
    weight_type = "rademacher",
    seed = 123L
  )

  expected <- contracts$monotonicity_surface
  expect_match(expected$source_anchor, "D19", fixed = TRUE)
  expect_equal(ti_result$att, expected$att, tolerance = 1e-12)
  expect_equal(pvalue_at_att, expected$pvalue_at_att, tolerance = 1e-12)
  expect_equal(pvalue_far_upper, expected$pvalue_far_upper, tolerance = 1e-12)
  expect_equal(pvalue_far_lower, expected$pvalue_far_lower, tolerance = 1e-12)
  expect_gt(pvalue_at_att, pvalue_far_upper)
  expect_gt(pvalue_at_att, pvalue_far_lower)
})

test_that("Task E9-06.13.11 E7 keeps the exact independent public test_inversion CI surface", {
  expect_true(
    exists("make_e906_public_boundary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.11 helper: make_e906_public_boundary_package_contracts"
  )
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_public_boundary_package_contracts", mode = "function") ||
    !exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_public_boundary_package_contracts()
  panel_data <- make_e906_public_lwdid_panel()
  ti_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = panel_data,
      y_transformed = "y",
      d = "treat",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      seed = 42L,
      grid_points = 11L,
      ci_tol = 0.01
    )
  )

  expected <- contracts$public_test_inversion
  expect_match(expected$source_anchor, "E7", fixed = TRUE)
  expect_identical(ti_result$ci_method, expected$ci_method)
  expect_identical(ti_result$actual_n_bootstrap, expected$actual_n_bootstrap)
  expect_equal(ti_result$att, expected$att, tolerance = 1e-12)
  expect_equal(ti_result$pvalue, expected$pvalue, tolerance = 1e-12)
  expect_equal(ti_result$ci_lower, expected$ci_lower, tolerance = 1e-12)
  expect_equal(ti_result$ci_upper, expected$ci_upper, tolerance = 1e-12)
  expect_equal(ti_result$original_se, expected$original_se, tolerance = 1e-12)
  expect_lte(ti_result$ci_lower, ti_result$att)
  expect_gte(ti_result$ci_upper, ti_result$att)
})

test_that("Task E9-06.13.11 F5 keeps the exact controls-null native same-path WCB surface", {
  expect_true(
    exists("make_e906_public_boundary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.11 helper: make_e906_public_boundary_package_contracts"
  )
  expect_true(
    exists("make_e906_restricted_model_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_restricted_model_panel"
  )

  if (!exists("make_e906_public_boundary_package_contracts", mode = "function") ||
    !exists("make_e906_restricted_model_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_public_boundary_package_contracts()
  restricted_model_panel <- make_e906_restricted_model_panel()
  with_controls_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = restricted_model_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      use_fwildclusterboot = FALSE,
      restricted_model = "with_controls"
    )
  )
  intercept_only_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = restricted_model_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      use_fwildclusterboot = FALSE,
      restricted_model = "intercept_only"
    )
  )

  expected <- contracts$controls_null_same_path
  expect_match(expected$source_anchor, "F5", fixed = TRUE)
  expect_identical(with_controls_result$actual_n_bootstrap, expected$actual_n_bootstrap)
  expect_identical(intercept_only_result$actual_n_bootstrap, expected$actual_n_bootstrap)
  expect_equal(with_controls_result$att, expected$att, tolerance = 1e-12)
  expect_equal(with_controls_result$pvalue, expected$pvalue, tolerance = 1e-12)
  expect_equal(with_controls_result$ci_lower, expected$ci_lower, tolerance = 1e-12)
  expect_equal(with_controls_result$ci_upper, expected$ci_upper, tolerance = 1e-12)
  expect_equal(intercept_only_result$att, expected$att, tolerance = 1e-12)
  expect_equal(intercept_only_result$pvalue, expected$pvalue, tolerance = 1e-12)
  expect_equal(intercept_only_result$ci_lower, expected$ci_lower, tolerance = 1e-12)
  expect_equal(intercept_only_result$ci_upper, expected$ci_upper, tolerance = 1e-12)
  expect_equal(with_controls_result$att, intercept_only_result$att, tolerance = 1e-12)
  expect_equal(with_controls_result$pvalue, intercept_only_result$pvalue, tolerance = 1e-12)
  expect_equal(with_controls_result$ci_lower, intercept_only_result$ci_lower, tolerance = 1e-12)
  expect_equal(with_controls_result$ci_upper, intercept_only_result$ci_upper, tolerance = 1e-12)
})

test_that("Task E9-06.13.11 F7 keeps the exact degenerate zero-SE test_inversion surface", {
  expect_true(
    exists("make_e906_public_boundary_package_contracts", mode = "function"),
    info = "missing Task E9-06.13.11 helper: make_e906_public_boundary_package_contracts"
  )
  expect_true(
    exists("make_e906_constant_outcome_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_constant_outcome_panel"
  )

  if (!exists("make_e906_public_boundary_package_contracts", mode = "function") ||
    !exists("make_e906_constant_outcome_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_public_boundary_package_contracts()
  constant_panel <- make_e906_constant_outcome_panel()
  ti_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = constant_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 49L,
      seed = 42L,
      grid_points = 5L
    )
  )

  expected <- contracts$degenerate_zero_se
  expect_match(expected$source_anchor, "F7", fixed = TRUE)
  expect_identical(ti_result$ci_method, expected$ci_method)
  expect_identical(ti_result$actual_n_bootstrap, expected$actual_n_bootstrap)
  expect_lte(abs(ti_result$att - expected$att), 1e-18)
  expect_lte(abs(ti_result$original_se - expected$original_se), 1e-18)
  expect_true(is.nan(ti_result$pvalue))
  expect_true(is.nan(ti_result$rejection_rate))
  expect_true(is.nan(ti_result$ci_lower))
  expect_true(is.nan(ti_result$ci_upper))
  expect_true(is.nan(expected$pvalue))
  expect_true(is.nan(expected$rejection_rate))
  expect_true(is.nan(expected$ci_lower))
  expect_true(is.nan(expected$ci_upper))
})

test_that("D16 direct test_inversion CI contains ATT on the shared-controls fixture", {
  expect_true(
    exists("make_e906_test_inversion_shared_fixture", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_test_inversion_shared_fixture"
  )

  if (!exists("make_e906_test_inversion_shared_fixture", mode = "function")) {
    return(invisible(NULL))
  }

  fixture <- make_e906_test_inversion_shared_fixture()
  ti_result <- suppressWarnings(
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

  expect_s3_class(ti_result, "lwdid_wcb_result")
  expect_identical(ti_result$ci_method, "test_inversion")
  expect_true(is.finite(ti_result$att))
  expect_true(is.finite(ti_result$ci_lower))
  expect_true(is.finite(ti_result$ci_upper))
  expect_lte(ti_result$ci_lower, ti_result$att)
  expect_gte(ti_result$ci_upper, ti_result$att)
})

test_that("D17 test_inversion CI width stays within 50% of percentile_t on the shared-controls fixture", {
  expect_true(
    exists("make_e906_test_inversion_shared_fixture", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_test_inversion_shared_fixture"
  )

  if (!exists("make_e906_test_inversion_shared_fixture", mode = "function")) {
    return(invisible(NULL))
  }

  fixture <- make_e906_test_inversion_shared_fixture()
  test_inversion <- suppressWarnings(
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
  percentile_t <- suppressWarnings(
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

  expect_true(is.finite(test_inversion_width))
  expect_true(is.finite(percentile_t_width))
  expect_gt(test_inversion_width, 0)
  expect_gt(percentile_t_width, 0)
  expect_gte(width_ratio, 0.5)
  expect_lte(width_ratio, 1.5)
  expect_lte(relative_width_gap, 0.5)
})

test_that("D19 test_inversion p-value is highest at ATT and declines at far nulls", {
  expect_true(
    exists("make_e906_test_inversion_shared_fixture", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_test_inversion_shared_fixture"
  )

  if (!exists("make_e906_test_inversion_shared_fixture", mode = "function")) {
    return(invisible(NULL))
  }

  fixture <- make_e906_test_inversion_shared_fixture()
  ti_result <- suppressWarnings(
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

  pvalue_at_att <- lwdid:::.compute_bootstrap_pvalue_at_null(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    null_value = ti_result$att,
    att_original = ti_result$att,
    se_original = ti_result$original_se,
    n_bootstrap = 999L,
    weight_type = "rademacher",
    seed = 123L
  )
  pvalue_far_upper <- lwdid:::.compute_bootstrap_pvalue_at_null(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    null_value = ti_result$att + 5 * ti_result$original_se,
    att_original = ti_result$att,
    se_original = ti_result$original_se,
    n_bootstrap = 999L,
    weight_type = "rademacher",
    seed = 123L
  )
  pvalue_far_lower <- lwdid:::.compute_bootstrap_pvalue_at_null(
    data = fixture,
    y_transformed = "Y",
    d = "D",
    cluster_var = "cluster",
    controls = "x1",
    null_value = ti_result$att - 5 * ti_result$original_se,
    att_original = ti_result$att,
    se_original = ti_result$original_se,
    n_bootstrap = 999L,
    weight_type = "rademacher",
    seed = 123L
  )

  expect_gt(pvalue_at_att, 0.5)
  expect_lt(pvalue_far_upper, pvalue_at_att)
  expect_lt(pvalue_far_lower, pvalue_at_att)
})

test_that("E7 independent public test_inversion call reports a finite test_inversion CI", {
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  panel_data <- make_e906_public_lwdid_panel()
  ti_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = panel_data,
      y_transformed = "y",
      d = "treat",
      cluster_var = "cluster",
      n_bootstrap = 99L,
      weight_type = "rademacher",
      seed = 42L,
      grid_points = 11L,
      ci_tol = 0.01
    )
  )

  expect_s3_class(ti_result, "lwdid_wcb_result")
  expect_identical(ti_result$ci_method, "test_inversion")
  expect_true(is.finite(ti_result$att))
  expect_true(is.finite(ti_result$pvalue))
  expect_true(is.finite(ti_result$ci_lower))
  expect_true(is.finite(ti_result$ci_upper))
  expect_lte(ti_result$ci_lower, ti_result$att)
  expect_gte(ti_result$ci_upper, ti_result$att)
})

test_that("F5 controls=NULL keeps the native WCB surface aligned across restricted models", {
  expect_true(
    exists("make_e906_public_lwdid_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_public_lwdid_panel"
  )

  if (!exists("make_e906_public_lwdid_panel", mode = "function")) {
    return(invisible(NULL))
  }

  panel_data <- make_e906_public_lwdid_panel()
  with_controls_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = panel_data,
      y_transformed = "y",
      d = "treat",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 99L,
      weight_type = "rademacher",
      seed = 42L,
      impose_null = TRUE,
      restricted_model = "with_controls",
      use_fwildclusterboot = FALSE
    )
  )
  intercept_only_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = panel_data,
      y_transformed = "y",
      d = "treat",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 99L,
      weight_type = "rademacher",
      seed = 42L,
      impose_null = TRUE,
      restricted_model = "intercept_only",
      use_fwildclusterboot = FALSE
    )
  )

  expect_s3_class(with_controls_result, "lwdid_wcb_result")
  expect_s3_class(intercept_only_result, "lwdid_wcb_result")
  expect_true(is.finite(with_controls_result$att))
  expect_true(is.finite(with_controls_result$pvalue))
  expect_equal(with_controls_result$att, intercept_only_result$att, tolerance = 1e-12)
  expect_equal(with_controls_result$pvalue, intercept_only_result$pvalue, tolerance = 1e-12)
  expect_equal(with_controls_result$ci_lower, intercept_only_result$ci_lower, tolerance = 1e-12)
  expect_equal(with_controls_result$ci_upper, intercept_only_result$ci_upper, tolerance = 1e-12)
})

test_that("F7 test_inversion returns a degenerate NaN CI when the constant-outcome fixture has zero usable SE", {
  expect_true(
    exists("make_e906_constant_outcome_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_constant_outcome_panel"
  )

  if (!exists("make_e906_constant_outcome_panel", mode = "function")) {
    return(invisible(NULL))
  }

  constant_panel <- make_e906_constant_outcome_panel()
  ti_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = constant_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 49L,
      seed = 42L,
      grid_points = 5L
    )
  )

  expect_s3_class(ti_result, "lwdid_wcb_result")
  expect_identical(ti_result$ci_method, "test_inversion")
  expect_true(is.finite(ti_result$att))
  expect_true(is.finite(ti_result$original_se))
  expect_lt(abs(ti_result$original_se), 1e-10)
  expect_true(is.nan(ti_result$pvalue))
  expect_true(is.nan(ti_result$rejection_rate))
  expect_true(is.nan(ti_result$ci_lower))
  expect_true(is.nan(ti_result$ci_upper))
})

test_that("D13 fwildclusterboot stays approximately aligned with the native controls-null package surface", {
  expect_true(
    exists("make_e906_test_inversion_panel", mode = "function"),
    info = "missing later-slice helper: make_e906_test_inversion_panel"
  )

  if (!exists("make_e906_test_inversion_panel", mode = "function")) {
    return(invisible(NULL))
  }

  skip_if_not_installed("fwildclusterboot")

  test_inversion_panel <- make_e906_test_inversion_panel()
  native_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = test_inversion_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      full_enumeration = FALSE,
      restricted_model = "intercept_only",
      use_fwildclusterboot = FALSE
    )
  )
  fwildclusterboot_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = test_inversion_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 99L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      full_enumeration = FALSE,
      restricted_model = "intercept_only",
      use_fwildclusterboot = TRUE
    )
  )

  expect_s3_class(native_result, "lwdid_wcb_result")
  expect_s3_class(fwildclusterboot_result, "lwdid_wcb_result")
  expect_true(is.finite(native_result$att))
  expect_true(is.finite(fwildclusterboot_result$att))
  expect_true(is.finite(native_result$pvalue))
  expect_true(is.finite(fwildclusterboot_result$pvalue))
  expect_equal(native_result$att, fwildclusterboot_result$att, tolerance = 1e-8)
  expect_equal(native_result$pvalue, fwildclusterboot_result$pvalue, tolerance = 0.05)
})

test_that("F5 controls=NULL keeps wild_cluster_bootstrap on the same native restricted-model path", {
  restricted_model_panel <- make_e906_restricted_model_panel()

  with_controls_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = restricted_model_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      use_fwildclusterboot = FALSE,
      restricted_model = "with_controls"
    )
  )
  intercept_only_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap(
      data = restricted_model_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 19L,
      weight_type = "rademacher",
      impose_null = TRUE,
      alpha = 0.05,
      seed = 42L,
      use_fwildclusterboot = FALSE,
      restricted_model = "intercept_only"
    )
  )

  expect_s3_class(with_controls_result, "lwdid_wcb_result")
  expect_s3_class(intercept_only_result, "lwdid_wcb_result")
  expect_true(is.finite(with_controls_result$att))
  expect_true(is.finite(with_controls_result$pvalue))
  expect_true(is.finite(with_controls_result$ci_lower))
  expect_true(is.finite(with_controls_result$ci_upper))
  expect_equal(with_controls_result$att, intercept_only_result$att, tolerance = 1e-12)
  expect_equal(with_controls_result$pvalue, intercept_only_result$pvalue, tolerance = 1e-12)
  expect_equal(with_controls_result$ci_lower, intercept_only_result$ci_lower, tolerance = 1e-12)
  expect_equal(with_controls_result$ci_upper, intercept_only_result$ci_upper, tolerance = 1e-12)
})

test_that("D16 test_inversion CI contains the ATT on the later native seasonal-WCB panel", {
  expect_true(
    exists("make_e906_test_inversion_panel", mode = "function"),
    info = "missing later-slice helper: make_e906_test_inversion_panel"
  )

  if (!exists("make_e906_test_inversion_panel", mode = "function")) {
    return(invisible(NULL))
  }

  test_inversion_panel <- make_e906_test_inversion_panel()
  result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = test_inversion_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 49L,
      seed = 42L,
      grid_points = 11L,
      ci_tol = 0.01
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_identical(result$ci_method, "test_inversion")
  expect_lte(result$ci_lower, result$att)
  expect_gte(result$ci_upper, result$att)
})

test_that("E7 standalone test_inversion advertises the native ci_method and finite interval", {
  expect_true(
    exists("make_e906_test_inversion_panel", mode = "function"),
    info = "missing later-slice helper: make_e906_test_inversion_panel"
  )

  if (!exists("make_e906_test_inversion_panel", mode = "function")) {
    return(invisible(NULL))
  }

  test_inversion_panel <- make_e906_test_inversion_panel()
  ti_result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = test_inversion_panel,
      y_transformed = "y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 49L,
      seed = 42L,
      grid_points = 11L,
      ci_tol = 0.01
    )
  )

  expect_s3_class(ti_result, "lwdid_wcb_result")
  expect_identical(ti_result$ci_method, "test_inversion")
  expect_true(is.finite(ti_result$pvalue))
  expect_true(is.finite(ti_result$ci_lower))
  expect_true(is.finite(ti_result$ci_upper))
  expect_lt(ti_result$ci_lower, ti_result$ci_upper)
})

test_that("F6 small-sample test_inversion keeps a finite ordered CI at N=12 and G=3", {
  small_ti_panel <- make_e906_small_ti_panel()
  result <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = small_ti_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      controls = NULL,
      n_bootstrap = 49L,
      seed = 42L,
      grid_points = 7L,
      ci_tol = 0.05
    )
  )

  expect_s3_class(result, "lwdid_wcb_result")
  expect_identical(result$ci_method, "test_inversion")
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
  expect_lt(result$ci_lower, result$ci_upper)
})

test_that("F10 test_inversion grid resolution converges to similar CI bounds", {
  expect_true(
    exists("make_e906_grid_resolution_panel", mode = "function"),
    info = "missing Task E9-06 helper: make_e906_grid_resolution_panel"
  )

  if (!exists("make_e906_grid_resolution_panel", mode = "function")) {
    return(invisible(NULL))
  }

  grid_panel <- make_e906_grid_resolution_panel()
  result_11 <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = grid_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 199L,
      seed = 42L,
      grid_points = 11L,
      ci_tol = 0.01
    )
  )
  result_21 <- suppressWarnings(
    lwdid::wild_cluster_bootstrap_test_inversion(
      data = grid_panel,
      y_transformed = "Y",
      d = "D",
      cluster_var = "cluster",
      n_bootstrap = 199L,
      seed = 42L,
      grid_points = 21L,
      ci_tol = 0.01
    )
  )

  lower_relative_gap <- abs(result_11$ci_lower - result_21$ci_lower) /
    max(abs(result_21$ci_lower), 1e-8)
  upper_relative_gap <- abs(result_11$ci_upper - result_21$ci_upper) /
    max(abs(result_21$ci_upper), 1e-8)

  expect_s3_class(result_11, "lwdid_wcb_result")
  expect_s3_class(result_21, "lwdid_wcb_result")
  expect_identical(result_11$ci_method, "test_inversion")
  expect_identical(result_21$ci_method, "test_inversion")
  expect_true(is.finite(result_11$ci_lower))
  expect_true(is.finite(result_11$ci_upper))
  expect_true(is.finite(result_21$ci_lower))
  expect_true(is.finite(result_21$ci_upper))
  expect_lte(lower_relative_gap, 0.05)
  expect_lte(upper_relative_gap, 0.05)
})

test_that("Layer 3 broader hardening demeanq smoking quarterly replay preserves the frozen realdata contract", {
  expect_true(
    exists("make_e906_layer3_realdata_contracts", mode = "function"),
    info = "missing owner helper: make_e906_layer3_realdata_contracts"
  )
  expect_true(
    exists("make_e906_smoking_quarterly_realdata_panel", mode = "function"),
    info = "missing owner helper: make_e906_smoking_quarterly_realdata_panel"
  )

  if (!exists("make_e906_layer3_realdata_contracts", mode = "function") ||
    !exists("make_e906_smoking_quarterly_realdata_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_layer3_realdata_contracts()
  panel_data <- make_e906_smoking_quarterly_realdata_panel()
  warnings_seen <- character(0)

  demeanq_result <- withCallingHandlers(
    lwdid::lwdid(
      data = panel_data,
      y = "lcigsale",
      d = "d",
      ivar = "state",
      tvar = c("year", "quarter"),
      post = "post",
      rolling = "demeanq",
      vce = NULL
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expected <- contracts$quarterly$cases$demeanq$expected
  expect_s3_class(demeanq_result, "lwdid_result")
  expect_equal(demeanq_result$att, expected$att, tolerance = 1e-12)
  expect_equal(demeanq_result$se_att, expected$se_att, tolerance = 1e-12)
  expect_equal(demeanq_result$ci_lower, expected$ci_lower, tolerance = 1e-12)
  expect_equal(demeanq_result$ci_upper, expected$ci_upper, tolerance = 1e-12)
  expect_equal(demeanq_result$pvalue, expected$pvalue, tolerance = 1e-12)
  expect_identical(demeanq_result$K, expected$K)
  expect_identical(demeanq_result$tpost1, expected$tpost1)
  expect_identical(nrow(demeanq_result$att_by_period), expected$rows)
  expect_identical(as.character(demeanq_result$att_by_period$period[[2L]]), expected$first_period)
  expect_identical(
    as.character(demeanq_result$att_by_period$period[[nrow(demeanq_result$att_by_period)]]),
    expected$last_period
  )
  expect_true(any(grepl(contracts$quarterly$cases$demeanq$warning_markers, warnings_seen, fixed = TRUE)))
})

test_that("Layer 3 broader hardening detrendq smoking quarterly replay preserves the frozen realdata contract", {
  expect_true(
    exists("make_e906_layer3_realdata_contracts", mode = "function"),
    info = "missing owner helper: make_e906_layer3_realdata_contracts"
  )
  expect_true(
    exists("make_e906_smoking_quarterly_realdata_panel", mode = "function"),
    info = "missing owner helper: make_e906_smoking_quarterly_realdata_panel"
  )

  if (!exists("make_e906_layer3_realdata_contracts", mode = "function") ||
    !exists("make_e906_smoking_quarterly_realdata_panel", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_layer3_realdata_contracts()
  panel_data <- make_e906_smoking_quarterly_realdata_panel()
  warnings_seen <- character(0)

  detrendq_result <- withCallingHandlers(
    lwdid::lwdid(
      data = panel_data,
      y = "lcigsale",
      d = "d",
      ivar = "state",
      tvar = c("year", "quarter"),
      post = "post",
      rolling = "detrendq",
      vce = NULL
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expected <- contracts$quarterly$cases$detrendq$expected
  expect_s3_class(detrendq_result, "lwdid_result")
  expect_equal(detrendq_result$att, expected$att, tolerance = 1e-12)
  expect_equal(detrendq_result$se_att, expected$se_att, tolerance = 1e-12)
  expect_equal(detrendq_result$ci_lower, expected$ci_lower, tolerance = 1e-12)
  expect_equal(detrendq_result$ci_upper, expected$ci_upper, tolerance = 1e-12)
  expect_equal(detrendq_result$pvalue, expected$pvalue, tolerance = 1e-12)
  expect_identical(detrendq_result$K, expected$K)
  expect_identical(detrendq_result$tpost1, expected$tpost1)
  expect_identical(nrow(detrendq_result$att_by_period), expected$rows)
  expect_identical(as.character(detrendq_result$att_by_period$period[[2L]]), expected$first_period)
  expect_identical(
    as.character(detrendq_result$att_by_period$period[[nrow(detrendq_result$att_by_period)]]),
    expected$last_period
  )
  expect_true(any(grepl(contracts$quarterly$cases$detrendq$warning_markers, warnings_seen, fixed = TRUE)))
})

test_that("Layer 3 broader hardening small-G public WCB replay preserves the frozen R boundary", {
  expect_true(
    exists("make_e906_layer3_realdata_contracts", mode = "function"),
    info = "missing owner helper: make_e906_layer3_realdata_contracts"
  )
  expect_true(
    exists("make_e906_smoking_small_g_public_wcb_fixture", mode = "function"),
    info = "missing owner helper: make_e906_smoking_small_g_public_wcb_fixture"
  )

  if (!exists("make_e906_layer3_realdata_contracts", mode = "function") ||
    !exists("make_e906_smoking_small_g_public_wcb_fixture", mode = "function")) {
    return(invisible(NULL))
  }

  contracts <- make_e906_layer3_realdata_contracts()
  fixture <- make_e906_smoking_small_g_public_wcb_fixture()

  for (seed_name in c("seed_123", "seed_999")) {
    expected <- contracts$small_g$seeds[[seed_name]]$r
    result <- suppressWarnings(
      lwdid::wild_cluster_bootstrap(
        data = fixture,
        y_transformed = "y_dot",
        d = "treated",
        cluster_var = "state",
        controls = NULL,
        n_bootstrap = 999L,
        weight_type = "rademacher",
        alpha = 0.05,
        seed = expected$seed,
        impose_null = TRUE,
        full_enumeration = NULL
      )
    )

    expect_s3_class(result, "lwdid_wcb_result")
    expect_equal(result$att, expected$att, tolerance = 1e-12)
    expect_equal(result$original_se, expected$original_se, tolerance = 1e-12)
    expect_equal(result$ci_lower, expected$ci_lower, tolerance = 1e-12)
    expect_equal(result$ci_upper, expected$ci_upper, tolerance = 1e-12)
    expect_equal(result$pvalue, expected$pvalue, tolerance = 1e-12)
    expect_equal(result$t_stat_original, expected$t_stat_original, tolerance = 1e-12)
    expect_identical(result$actual_n_bootstrap, expected$actual_n_bootstrap)
    expect_identical(isTRUE(result$full_enumeration), expected$full_enumeration)
    expect_identical(result$ci_method, expected$ci_method)
  }
})
