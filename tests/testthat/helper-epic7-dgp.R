# ===========================================================================
# helper-epic7-dgp.R
# DGP helper functions for Epic 7 comprehensive test suite (Story E7-07)
# Automatically loaded by testthat before test execution.
# ===========================================================================

#' Generate Common Timing panel data
#' @param N Total units (default 100)
#' @param T_total Total periods (default 8)
#' @param tpost1 First post-treatment period (default 5, i.e. S=5)
#' @param treat_frac Treatment fraction (default 0.4)
#' @param tau True ATT (default 2.0)
#' @param pt_violation PT violation magnitude (default 0)
#' @param with_controls Include control variables (default FALSE)
#' @param n_controls Number of controls (default 2)
#' @param seed Random seed
#' @return data.table with id/time/y/d columns
generate_ct_panel <- function(N = 100L, T_total = 8L, tpost1 = 5L,
                              treat_frac = 0.4, tau = 2.0,
                              pt_violation = 0, with_controls = FALSE,
                              n_controls = 2L, seed = 42L) {
  set.seed(seed)
  n_treat <- as.integer(round(N * treat_frac))
  n_ctrl <- N - n_treat
  periods <- seq_len(T_total)

  dt <- data.table::CJ(id = seq_len(N), time = periods)
  dt[, d := as.integer(id <= n_treat)]

  # Unit FE ~ U(5,15), common time slope ~ U(0.5,1.5)
  alpha_i <- stats::runif(N, 5, 15)
  beta <- stats::runif(1, 0.5, 1.5)

  dt[, y := alpha_i[id] + beta * time +
       tau * d * as.integer(time >= tpost1) +
       stats::rnorm(.N, sd = 0.5)]

  # PT violation: treated units get extra linear trend
  if (pt_violation > 0) {
    dt[d == 1L, y := y + pt_violation * time]
  }

  # Control variables
  if (with_controls) {
    for (k in seq_len(n_controls)) {
      col_name <- paste0("x", k)
      x_vals <- stats::rnorm(N)
      dt[, (col_name) := x_vals[id]]
    }
    # Controls affect Y: 0.3*x1 - 0.2*x2
    if (n_controls >= 1L) dt[, y := y + 0.3 * x1]
    if (n_controls >= 2L) dt[, y := y - 0.2 * x2]
  }
  dt
}


#' Generate Staggered panel data
#' @param n_per_cohort Units per cohort (default 5)
#' @param cohorts Cohort vector (default c(4,6,8))
#' @param n_never_treated Never-treated units (default 5)
#' @param T_total Total periods (default 10)
#' @param tau_base Base treatment effect (default 2.0)
#' @param tau_dynamic Dynamic effect slope (default 0.5)
#' @param pt_violation PT violation magnitude (default 0)
#' @param with_controls Include controls (default FALSE)
#' @param seed Random seed
#' @return data.table with id/time/y/gvar columns
generate_staggered_panel <- function(n_per_cohort = 5L,
                                     cohorts = c(4L, 6L, 8L),
                                     n_never_treated = 5L,
                                     T_total = 10L,
                                     tau_base = 2.0,
                                     tau_dynamic = 0.5,
                                     pt_violation = 0,
                                     with_controls = FALSE,
                                     seed = 42L) {
  set.seed(seed)
  n_treat <- length(cohorts) * n_per_cohort
  n_total <- n_treat + n_never_treated
  periods <- seq_len(T_total)

  # Cohort assignments: treated cohorts then never-treated (Inf)
  gvar_vals <- c(
    rep(cohorts, each = n_per_cohort),
    rep(Inf, n_never_treated)
  )

  dt <- data.table::CJ(id = seq_len(n_total), time = periods)
  dt[, gvar := gvar_vals[id]]

  # Unit FE ~ U(5,15), unit-specific time slope ~ U(0.5,1.5)
  alpha_i <- stats::runif(n_total, 5, 15)
  beta_i <- stats::runif(n_total, 0.5, 1.5)

  dt[, y := alpha_i[id] + beta_i[id] * time + stats::rnorm(.N, sd = 0.5)]

  # Treatment effect: tau(t,g) = tau_base + tau_dynamic*(t-g)

  dt[gvar < Inf & time >= gvar,
     y := y + tau_base + tau_dynamic * (time - gvar)]

  # PT violation: treated units get extra trend in pre-treatment
  if (pt_violation > 0) {
    dt[gvar < Inf & time < gvar,
       y := y + pt_violation * time]
  }

  # Control variables
  if (with_controls) {
    x_vals <- stats::rnorm(n_total)
    dt[, x1 := x_vals[id]]
    dt[, y := y + 0.3 * x1]
  }
  dt
}

#' Generate minimal panel (N=3, T=3)
#' @param seed Random seed
#' @return data.table
generate_minimal_panel <- function(seed = 42L) {
  set.seed(seed)
  data.table::data.table(
    id = rep(1:3, each = 3L),
    time = rep(1:3, 3L),
    y = stats::rnorm(9),
    d = rep(c(1L, 0L, 0L), each = 3L),
    gvar = rep(c(2L, Inf, Inf), each = 3L)
  )
}

#' Generate constant outcome panel (Y=5.0)
#' @param N Total units (default 20)
#' @param T_total Total periods (default 6)
#' @param seed Random seed
#' @return data.table
generate_constant_panel <- function(N = 20L, T_total = 6L, seed = 42L) {
  set.seed(seed)
  n_treat <- N %/% 2L
  dt <- data.table::CJ(id = seq_len(N), time = seq_len(T_total))
  dt[, d := as.integer(id <= n_treat)]
  dt[, y := 5.0]
  dt
}

#' Generate extreme value panel (Y ~ 1e10)
#' @param N Total units (default 20)
#' @param T_total Total periods (default 6)
#' @param seed Random seed
#' @return data.table
generate_extreme_panel <- function(N = 20L, T_total = 6L, seed = 42L) {
  set.seed(seed)
  n_treat <- N %/% 2L
  dt <- data.table::CJ(id = seq_len(N), time = seq_len(T_total))
  dt[, d := as.integer(id <= n_treat)]
  dt[, y := 1e10 + stats::rnorm(.N, sd = 1)]
  dt
}

#' Generate staggered panel with no never-treated units
#' @param seed Random seed
#' @return data.table
generate_no_nt_panel <- function(seed = 42L) {
  set.seed(seed)
  dt <- data.table::CJ(id = seq_len(12), time = seq_len(8))
  dt[, gvar := rep(c(3L, 5L, 7L), each = 4L)[id]]
  alpha_i <- stats::rnorm(12, sd = 2)
  dt[, y := alpha_i[id] + 0.5 * time + stats::rnorm(.N, sd = 0.3)]
  dt[time >= gvar, y := y + 2.0]
  dt
}
