# Extracted from test-vce-unit.R:271

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
y <- c(3, 3, 3, 1, 1, 1)
d <- c(1, 1, 1, 0, 0, 0)
fit <- lm(y ~ d)
expect_warning(
    compute_vce(fit, vce = NULL),
    class = "lwdid_numerical"
  )
