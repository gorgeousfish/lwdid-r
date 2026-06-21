# Extracted from test-vce-unit.R:690

# prequel ----------------------------------------------------------------------
generate_vce_test_data <- function(n = 30, n_treated = 15, k = 2,
                                   seed = 42) {
  set.seed(seed)
  d <- c(rep(1L, n_treated), rep(0L, n - n_treated))
  x <- if (k > 0) matrix(rnorm(n * k), n, k) else NULL
  y <- 1 + 2 * d + (if (k > 0) rowSums(x) else 0) + rnorm(n)
  list(y = y, d = d, x = x)
}
generate_clustered_data <- function(g = 10, obs_per_cluster = 6,
                                    seed = 42) {
  set.seed(seed)
  if (length(obs_per_cluster) == 1) {
    obs_per_cluster <- rep(obs_per_cluster, g)
  }
  n <- sum(obs_per_cluster)
  cluster <- rep(seq_len(g), times = obs_per_cluster)
  # Cluster-level random effect
  cluster_effect <- rnorm(g, sd = 2)[cluster]
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + cluster_effect + rnorm(n)
  list(y = y, d = d, x = NULL, cluster = cluster, n = n, g = g)
}

# test -------------------------------------------------------------------------
set.seed(42)
n <- 20
d <- c(rep(1, 10), rep(0, 10))
x <- rnorm(n)
y <- 1 + 2 * d + 0.5 * x + rnorm(n)
fit <- lm(y ~ d + x)
V_hc3 <- sandwich::vcovHC(fit, type = "HC3")
beta_full <- coef(fit)
beta_jack <- matrix(NA, nrow = n, ncol = length(beta_full))
for (i in seq_len(n)) {
    fit_i <- lm(y ~ d + x, subset = -i)
    beta_jack[i, ] <- coef(fit_i)
  }
V_jack <- ((n - 1) / n) *
    crossprod(sweep(beta_jack, 2, beta_full, "-"))
V_expected <- ((n - 1) / n) * V_hc3
expect_equal(V_jack, V_expected, tolerance = 1e-10,
               label = "jackknife variance = (N-1)/N * HC3 variance")
