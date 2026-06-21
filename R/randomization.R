#' @title Randomization Inference for ATT
#'
#' @description Performs randomization inference (RI) to compute
#'   p-values for the average treatment effect on the treated (ATT).
#'   Supports both bootstrap and permutation methods, with optional
#'   covariate adjustment following Liang & Wang (2025) equations
#'   3.3 and 3.4.
#'
#' @param y_trans Numeric vector of transformed outcomes (length N).
#' @param d Integer vector of treatment indicators (0 or 1, length N).
#' @param x Optional numeric matrix of covariates (N rows). When
#'   provided, the ATT is estimated via OLS with interactions
#'   centered at the treatment group mean (lw2025 eq 3.3).
#' @param reps Integer number of RI replications (default 1000).
#' @param seed Optional integer seed for reproducibility.
#' @param method Character; either \code{"bootstrap"} (default) or
#'   \code{"permutation"}.
#' @param alpha Numeric significance level (default 0.05).
#'
#' @return A list with 9 fields:
#'   \describe{
#'     \item{ri_pvalue}{Two-sided RI p-value.}
#'     \item{ri_distribution}{Numeric vector of RI test statistics.}
#'     \item{obs_att}{Observed ATT estimate.}
#'     \item{reps}{Number of replications requested.}
#'     \item{n_valid}{Number of valid (non-failed) replications.}
#'     \item{n_failed}{Number of failed replications.}
#'     \item{failure_rate}{Proportion of failed replications.}
#'     \item{method}{Character string of the method used.}
#'     \item{seed}{Seed used, or NULL.}
#'   }
#'
#' @keywords internal
randomization_inference <- function(y_trans, d, x = NULL,
                                    reps = 1000L, seed = NULL,
                                    method = c("bootstrap",
                                               "permutation"),
                                    alpha = 0.05,
                                    parallel = FALSE) {

  # -- Parameter matching -----------------------------------------------------
  method <- match.arg(method)

  # -- Seed handling ----------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  # -- Input validation -------------------------------------------------------
  if (length(y_trans) != length(d)) {
    stop_lwdid(
      "y_trans and d must have the same length",
      class = "lwdid_invalid_param"
    )
  }

  if (!all(d %in% c(0L, 1L))) {
    stop_lwdid(
      "d must contain only 0 and 1 values",
      class = "lwdid_invalid_param"
    )
  }

  if (anyNA(y_trans)) {
    stop_lwdid(
      "y_trans must not contain NA values",
      class = "lwdid_invalid_param"
    )
  }

  if (!is.null(x)) {
    if (!is.matrix(x)) x <- as.matrix(x)
    if (nrow(x) != length(y_trans)) {
      stop_lwdid(
        "x must have the same number of rows as length(y_trans)",
        class = "lwdid_invalid_param"
      )
    }
    if (anyNA(x)) {
      stop_lwdid(
        "x must not contain NA values",
        class = "lwdid_invalid_param"
      )
    }
  }

  if (reps < 1L) {
    stop_lwdid(
      "reps must be a positive integer (>= 1)",
      class = "lwdid_invalid_param"
    )
  }

  # -- Basic quantities -------------------------------------------------------
  n  <- length(y_trans)
  n1 <- sum(d == 1L)

  # -- Sample size check ------------------------------------------------------
  if (n < 3L) {
    stop_lwdid(
      paste0(
        "Sample size too small for randomization inference: ",
        "N must be >= 3"
      ),
      class = "lwdid_insufficient_sample"
    )
  }

  # -- Observed ATT calculation -----------------------------------------------
  if (is.null(x)) {
    # Without controls (lw2025 eq 3.4)
    obs_att <- mean(y_trans[d == 1L]) - mean(y_trans[d == 0L])
  } else {
    # With controls (lw2025 eq 3.3, FATAL-003 protection)
    # Center covariates at treatment group mean (FIXED, outside loop)
    x_mean_treat <- colMeans(x[d == 1L, , drop = FALSE])
    x_c <- sweep(x, 2, x_mean_treat)
    interactions <- d * x_c
    X_mat <- cbind(1, d, x_c, interactions)
    qr_obs <- qr(X_mat)
    coefs <- qr.coef(qr_obs, y_trans)
    obs_att <- unname(coefs[2L])
  }

  # --- RI loop (Task E7-01.2) ---
  # Closure for a single RI replication
  .ri_one_rep <- function(rep_idx) {
    if (method == "permutation") {
      d_perm <- sample(d)
    } else {
      d_perm <- d[sample.int(n, n, replace = TRUE)]
      n1_perm <- sum(d_perm == 1L)
      if (n1_perm == 0L || n1_perm == n) return(NA_real_)
    }
    if (is.null(x)) {
      mean(y_trans[d_perm == 1L]) - mean(y_trans[d_perm == 0L])
    } else {
      inter_perm <- d_perm * x_c
      X_perm <- cbind(1, d_perm, x_c, inter_perm)
      qr_perm <- qr(X_perm)
      if (qr_perm$rank < ncol(X_perm)) return(NA_real_)
      qr.coef(qr_perm, y_trans)[2L]
    }
  }

  # Convergence-monitored RI loop: track failures in real time
  ri_early_stop <- FALSE

  if (isTRUE(parallel) && .can_use_parallel()) {
    ri_stats <- unlist(
      .parallel_lapply(seq_len(reps), .ri_one_rep),
      use.names = FALSE
    )
  } else {
    # Sequential path with convergence monitoring
    ri_stats <- rep(NA_real_, reps)
    n_failed_running <- 0L

    for (b in seq_len(reps)) {
      ri_stats[b] <- .ri_one_rep(b)

      if (is.na(ri_stats[b])) {
        n_failed_running <- n_failed_running + 1L
      }

      # Convergence check: after 50 iterations, if >20% failed, warn and stop
      if (b >= 50L && n_failed_running / b > 0.2) {
        warn_lwdid(
          sprintf(
            paste0("Randomization inference: %.1f%% of permutations failed ",
                   "(%d/%d). Results may be unreliable. ",
                   "Consider checking data structure or reducing reps."),
            n_failed_running / b * 100, n_failed_running, b
          ),
          class = "lwdid_ri_convergence"
        )
        ri_stats <- ri_stats[seq_len(b)]
        ri_early_stop <- TRUE
        break
      }
    }
  }

  # --- Post-processing (Task E7-01.3) ---

  # 1. Remove failed (NA) replications
  ri_valid <- ri_stats[!is.na(ri_stats)]

  # 2. Compute failure statistics
  n_total      <- length(ri_stats)
  n_valid      <- length(ri_valid)
  n_failed     <- n_total - n_valid
  failure_rate <- n_failed / max(n_total, 1L)

  # 3. Final high-failure-rate warning (post-loop, complements early-stop)
  if (!ri_early_stop && failure_rate > 0.1) {
    warn_lwdid(
      sprintf(
        paste0("Randomization inference completed with %.1f%% failure rate. ",
               "Effective replications: %d/%d"),
        failure_rate * 100, n_valid, reps
      ),
      class = "lwdid_ri_convergence"
    )
  } else if (!ri_early_stop && failure_rate > 0.05) {
    warn_lwdid(
      sprintf(
        "RI: %d/%d replications failed (%.1f%%); results may be unreliable",
        n_failed, reps, failure_rate * 100
      ),
      class = "lwdid_numerical"
    )
  }

  # 4. Minimum valid count threshold (method-dependent)
  if (method == "bootstrap") {
    min_valid <- max(100L, as.integer(0.1 * reps))
  } else {
    min_valid <- max(10L, as.integer(0.1 * reps))
  }

  # 5. Error if insufficient valid replications
  if (n_valid < min_valid) {
    stop_lwdid(
      sprintf(
        "RI: only %d valid replications (minimum %d required)",
        n_valid, min_valid
      ),
      class = "lwdid_numerical"
    )
  }

  # 6. Two-sided p-value as the simple randomization tail proportion.
  ri_pvalue <- sum(abs(ri_valid) >= abs(obs_att)) / n_valid

  # 7. Return result list (9 fields)
  list(
    ri_pvalue       = ri_pvalue,
    ri_distribution = ri_valid,
    obs_att         = obs_att,
    reps            = reps,
    n_valid         = n_valid,
    n_failed        = n_failed,
    failure_rate    = failure_rate,
    method          = method,
    seed            = seed
  )
}

.common_timing_randomization_inference <- function(
    data, y, d, controls = NULL, ps_controls = NULL,
    estimator = c("ipw", "ipwra", "psm"),
    reps = 1000L, seed = NULL,
    method = c("bootstrap", "permutation"),
    trim_threshold = 0.01,
    trim_method = "clip",
    n_neighbors = 1L,
    caliper = NULL,
    caliper_scale = "sd",
    with_replacement = TRUE,
    match_order = "data",
    se_method = NULL,
    boot_reps = 200L,
    return_diagnostics = FALSE,
    parallel = FALSE) {

  estimator <- match.arg(estimator)
  method <- match.arg(method)

  if (!is.null(seed)) set.seed(seed)

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  if (!y %in% names(data) || !d %in% names(data)) {
    stop_lwdid(
      "Common-timing RI requires transformed outcome and treatment columns in data.",
      class = "lwdid_invalid_param"
    )
  }
  if (anyNA(data[[y]]) || anyNA(data[[d]])) {
    stop_lwdid(
      "Common-timing RI data must not contain NA transformed outcomes or treatment values.",
      class = "lwdid_invalid_param"
    )
  }
  d_vec <- as.integer(data[[d]])
  if (!all(d_vec %in% c(0L, 1L))) {
    stop_lwdid(
      "Common-timing RI treatment column must contain only 0 and 1 values.",
      class = "lwdid_invalid_param"
    )
  }
  reps <- as.integer(reps)
  if (reps < 1L) {
    stop_lwdid(
      "reps must be a positive integer (>= 1)",
      class = "lwdid_invalid_param"
    )
  }

  n <- length(d_vec)
  if (n < 3L) {
    stop_lwdid(
      paste0(
        "Sample size too small for common-timing randomization inference: ",
        "N must be >= 3"
      ),
      class = "lwdid_insufficient_sample"
    )
  }

  estimate_att <- function(d_candidate) {
    est_data <- data
    est_data[[d]] <- as.integer(d_candidate)
    res <- dispatch_estimator(
      data = est_data,
      y = y,
      d = d,
      controls = controls,
      ps_controls = ps_controls,
      estimator = estimator,
      vce = NULL,
      cluster_var = NULL,
      alpha = 0.05,
      trim_threshold = trim_threshold,
      trim_method = trim_method,
      n_neighbors = n_neighbors,
      caliper = caliper,
      caliper_scale = caliper_scale,
      with_replacement = with_replacement,
      match_order = match_order,
      se_method = se_method,
      boot_reps = boot_reps,
      seed = seed,
      return_diagnostics = return_diagnostics
    )
    as.numeric(res$att)
  }

  obs_att <- estimate_att(d_vec)

  .ri_one_rep <- function(rep_idx) {
    d_perm <- if (method == "permutation") {
      sample(d_vec)
    } else {
      d_vec[sample.int(n, n, replace = TRUE)]
    }
    n1_perm <- sum(d_perm == 1L)
    if (n1_perm == 0L || n1_perm == n) {
      return(NA_real_)
    }
    tryCatch(
      estimate_att(d_perm),
      error = function(e) NA_real_
    )
  }

  # Convergence-monitored RI loop
  ri_early_stop <- FALSE

  if (isTRUE(parallel) && .can_use_parallel()) {
    ri_stats <- unlist(
      .parallel_lapply(seq_len(reps), .ri_one_rep),
      use.names = FALSE
    )
  } else {
    # Sequential path with convergence monitoring
    ri_stats <- rep(NA_real_, reps)
    n_failed_running <- 0L

    for (b in seq_len(reps)) {
      ri_stats[b] <- .ri_one_rep(b)

      if (is.na(ri_stats[b])) {
        n_failed_running <- n_failed_running + 1L
      }

      # Convergence check: after 50 iterations, if >20% failed, warn and stop
      if (b >= 50L && n_failed_running / b > 0.2) {
        warn_lwdid(
          sprintf(
            paste0("Randomization inference: %.1f%% of permutations failed ",
                   "(%d/%d). Results may be unreliable. ",
                   "Consider checking data structure or reducing reps."),
            n_failed_running / b * 100, n_failed_running, b
          ),
          class = "lwdid_ri_convergence"
        )
        ri_stats <- ri_stats[seq_len(b)]
        ri_early_stop <- TRUE
        break
      }
    }
  }

  ri_valid <- ri_stats[!is.na(ri_stats)]
  n_total <- length(ri_stats)
  n_valid <- length(ri_valid)
  n_failed <- n_total - n_valid
  failure_rate <- n_failed / max(n_total, 1L)

  # Final high-failure-rate warning (post-loop, complements early-stop)
  if (!ri_early_stop && failure_rate > 0.1) {
    warn_lwdid(
      sprintf(
        paste0("Randomization inference completed with %.1f%% failure rate. ",
               "Effective replications: %d/%d"),
        failure_rate * 100, n_valid, reps
      ),
      class = "lwdid_ri_convergence"
    )
  } else if (!ri_early_stop && failure_rate > 0.05) {
    warn_lwdid(
      sprintf(
        "RI: %d/%d replications failed (%.1f%%); results may be unreliable",
        n_failed, reps, failure_rate * 100
      ),
      class = "lwdid_numerical"
    )
  }

  min_valid <- if (method == "bootstrap") {
    max(100L, as.integer(0.1 * reps))
  } else {
    max(10L, as.integer(0.1 * reps))
  }
  if (n_valid < min_valid) {
    stop_lwdid(
      sprintf(
        "RI: only %d valid replications (minimum %d required)",
        n_valid, min_valid
      ),
      class = "lwdid_numerical"
    )
  }

  list(
    ri_pvalue = sum(abs(ri_valid) >= abs(obs_att)) / n_valid,
    ri_distribution = ri_valid,
    obs_att = obs_att,
    reps = reps,
    n_valid = n_valid,
    n_failed = n_failed,
    failure_rate = failure_rate,
    method = method,
    seed = seed
  )
}


#' @title Staggered Randomization Inference
#'
#' @description Performs randomization inference for staggered DiD
#'   by permuting unit-level cohort assignments. Each permutation
#'   re-runs the full pipeline (transform -> estimate -> aggregate)
#'   to build a reference distribution under the sharp null
#'   \eqn{H_0: \tau_i = 0, \forall i} (lw2026 Section 2).
#'
#' @param data data.table, original panel data.
#' @param y Character, outcome variable name.
#' @param ivar Character, unit identifier variable name.
#' @param tvar Character, time variable name.
#' @param gvar Character, cohort variable name.
#' @param observed_att Numeric, observed ATT from the caller.
#' @param target Character, inference target:
#'   \code{"overall"}, \code{"cohort"}, or \code{"cohort_time"}.
#' @param target_cohort Integer or NULL, target cohort
#'   (required when \code{target = "cohort"} or
#'   \code{"cohort_time"}).
#' @param target_period Integer or NULL, target period
#'   (required when \code{target = "cohort_time"}).
#' @param rolling Character, panel transformation method
#'   (\code{"demean"} or \code{"detrend"}).
#' @param estimator Character, estimator used when recomputing each
#'   permuted sample. One of \code{"ra"}, \code{"ipw"}, \code{"ipwra"},
#'   or \code{"psm"}.
#' @param controls Character vector or NULL, outcome-regression control
#'   names forwarded to estimation functions.
#' @param ps_controls Character vector or NULL, propensity-score control
#'   names for IPW/IPWRA/PSM permutation refits.
#' @param vce Character or NULL, variance estimation method.
#' @param cluster_var Character or NULL, clustering variable.
#' @param reps Positive integer, number of RI replications
#'   (default 1000).
#' @param seed Integer or NULL, random seed for reproducibility.
#' @param method Character, \code{"permutation"} (default) or
#'   \code{"bootstrap"}.
#' @param trim_threshold,trim_method,n_neighbors,caliper,caliper_scale,with_replacement,match_order,se_method,boot_reps
#'   Estimator-specific tuning parameters forwarded when
#'   \code{estimator} is IPW, IPWRA, or PSM.
#'
#' @return A list with 11 fields:
#'   \describe{
#'     \item{ri_pvalue}{Two-sided RI p-value.}
#'     \item{ri_distribution}{Numeric vector of valid RI stats.}
#'     \item{obs_att}{Observed ATT passed in.}
#'     \item{target}{Inference target used.}
#'     \item{target_cohort}{Target cohort, or NULL.}
#'     \item{target_period}{Target period, or NULL.}
#'     \item{n_valid}{Number of valid replications.}
#'     \item{n_failed}{Number of failed replications.}
#'     \item{failure_rate}{Proportion of failed replications.}
#'     \item{reps}{Number of replications requested.}
#'     \item{method}{Method used.}
#'   }
#'
#' @keywords internal
ri_staggered <- function(data, y, ivar, tvar, gvar,
                         observed_att,
                         target = c("overall", "cohort",
                                    "cohort_time"),
                         target_cohort = NULL,
                         target_period = NULL,
                         rolling = "demean",
                         estimator = "ra",
                         controls = NULL,
                         ps_controls = NULL,
                         vce = NULL, cluster_var = NULL,
                         reps = 1000L, seed = NULL,
                         method = c("permutation",
                                    "bootstrap"),
                         trim_threshold = 0.01,
                         trim_method = "clip",
                         n_neighbors = 1L,
                         caliper = NULL,
                         caliper_scale = "sd",
                         with_replacement = TRUE,
                         match_order = "data",
                         se_method = NULL,
                         boot_reps = 200L,
                         parallel = FALSE) {

  # -- Phase 1: Input validation & initialization --

  target <- match.arg(target)
  method <- match.arg(method)
  estimator <- match.arg(estimator, c("ra", "ipw", "ipwra", "psm"))

  if (!is.null(seed)) set.seed(seed)

  if (reps < 1L) {
    stop_lwdid(
      "reps must be a positive integer (>= 1)",
      class = "lwdid_invalid_param"
    )
  }

  if (target == "cohort" && is.null(target_cohort)) {
    stop_lwdid(
      "target_cohort is required when target = 'cohort'",
      class = "lwdid_invalid_param"
    )
  }

  if (target == "cohort_time" &&
      (is.null(target_cohort) || is.null(target_period))) {
    stop_lwdid(
      paste0("target_cohort and target_period are both ",
             "required when target = 'cohort_time'"),
      class = "lwdid_invalid_param"
    )
  }

  # -- Extract unit-level cohort assignments --
  unit_gvar_dt <- get_unit_level_gvar(
    data.table::as.data.table(data), gvar, ivar
  )
  gvar_values <- unit_gvar_dt[[gvar]]
  unit_ids    <- unit_gvar_dt[[ivar]]
  n_units     <- length(unit_ids)

  if (n_units < 4L) {
    stop_lwdid(
      sprintf(
        paste0("Staggered RI requires at least 4 units, ",
               "but only %d found"),
        n_units
      ),
      class = "lwdid_insufficient_sample"
    )
  }

  # For overall target, require never-treated units
  if (target == "overall") {
    n_nt <- sum(is_never_treated(gvar_values))
    if (n_nt == 0L) {
      stop_lwdid(
        paste0("Overall RI requires never-treated units, ",
               "but none found in data"),
        class = "lwdid_invalid_param"
      )
    }
  }

  t_max <- max(data[[tvar]], na.rm = TRUE)

  # -- Phase 2: RI loop --

  # Closure for a single staggered RI replication
  .ri_stag_one_rep <- function(b) {
    # Step 1: Permute unit-level cohort assignments
    if (method == "permutation") {
      perm_gvar <- gvar_values[sample.int(n_units)]
    } else {
      perm_gvar <- gvar_values[
        sample.int(n_units, n_units, replace = TRUE)
      ]
    }

    # Step 2: Build permuted data (replace gvar column)
    gvar_map  <- stats::setNames(perm_gvar, unit_ids)
    data_perm <- data.table::copy(data)
    data_perm[[gvar]] <- unname(
      gvar_map[as.character(data_perm[[ivar]])]
    )

    # Step 3: Re-run full pipeline in tryCatch
    tryCatch({
      perm_unit_gvar <- get_unit_level_gvar(
        data_perm, gvar, ivar
      )
      perm_gvar_vals <- perm_unit_gvar[[gvar]]
      perm_cohorts <- sort(unique(
        perm_gvar_vals[!is_never_treated(perm_gvar_vals)]
      ))
      perm_cohorts <- as.integer(perm_cohorts)

      if (length(perm_cohorts) == 0L) return(NA_real_)

      pre_stats <- suppressWarnings(
        precompute_transforms(
          data_perm, y, ivar, tvar,
          perm_cohorts, rolling
        )
      )
      valid_cohorts <- sort(as.integer(names(pre_stats)))

      if (length(valid_cohorts) == 0L) return(NA_real_)

      if (target == "overall") {
        agg <- suppressWarnings(
          if (identical(estimator, "ra")) {
            aggregate_to_overall(
              data_perm, y, ivar, tvar, gvar,
              valid_cohorts, t_max, pre_stats, rolling,
              vce = vce, cluster_var = cluster_var,
              controls = controls
            )
          } else {
            ct_effects <- estimate_staggered_effects(
              data_perm, y, ivar, tvar, gvar,
              rolling, "never_treated", controls,
              vce, cluster_var, alpha = 0.05,
              pre_stats,
              estimator = estimator,
              ps_controls = ps_controls,
              trim_threshold = trim_threshold,
              trim_method = trim_method,
              n_neighbors = n_neighbors,
              caliper = caliper,
              caliper_scale = caliper_scale,
              with_replacement = with_replacement,
              match_order = match_order,
              se_method = se_method,
              boot_reps = boot_reps
            )
            cohort_effects <- .aggregate_gr_to_cohort(
              ct_effects, valid_cohorts, alpha = 0.05
            )
            cohort_units <- get_unit_level_gvar(
              data.table::as.data.table(data_perm), gvar, ivar
            )
            cohort_sizes <- cohort_units[
              !is_never_treated(cohort_units[[gvar]]),
              ,
              drop = FALSE
            ]
            cohort_sizes <- data.table::as.data.table(cohort_sizes)[
              ,
              .N,
              by = gvar
            ]
            cohort_size_vec <- stats::setNames(
              cohort_sizes$N,
              as.character(cohort_sizes[[gvar]])
            )
            .aggregate_non_ra_cohorts_to_overall(
              cohort_effects = cohort_effects,
              cohort_sizes = cohort_size_vec,
              alpha = 0.05
            )
          }
        )
        return(agg$att)

      } else if (target == "cohort") {
        if (!(target_cohort %in% valid_cohorts)) return(NA_real_)
        cohort_res <- suppressWarnings(
          if (identical(estimator, "ra")) {
            aggregate_to_cohort(
              data_perm, y, ivar, tvar, gvar,
              target_cohort, t_max, pre_stats, rolling,
              vce = vce, cluster_var = cluster_var,
              controls = controls
            )
          } else {
            ct_effects <- estimate_staggered_effects(
              data_perm, y, ivar, tvar, gvar,
              rolling, "never_treated", controls,
              vce, cluster_var, alpha = 0.05,
              pre_stats,
              estimator = estimator,
              ps_controls = ps_controls,
              trim_threshold = trim_threshold,
              trim_method = trim_method,
              n_neighbors = n_neighbors,
              caliper = caliper,
              caliper_scale = caliper_scale,
              with_replacement = with_replacement,
              match_order = match_order,
              se_method = se_method,
              boot_reps = boot_reps
            )
            .aggregate_gr_to_cohort(
              ct_effects, target_cohort, alpha = 0.05
            )
          }
        )
        return(cohort_res[[1L]]$att)

      } else {
        ct_effects <- suppressWarnings(
          estimate_staggered_effects(
            data_perm, y, ivar, tvar, gvar,
            rolling, "never_treated", controls,
            vce, cluster_var, alpha = 0.05,
            pre_stats,
            estimator = estimator,
            ps_controls = ps_controls,
            trim_threshold = trim_threshold,
            trim_method = trim_method,
            n_neighbors = n_neighbors,
            caliper = caliper,
            caliper_scale = caliper_scale,
            with_replacement = with_replacement,
            match_order = match_order,
            se_method = se_method,
            boot_reps = boot_reps
          )
        )
        match_row <- ct_effects[
          ct_effects$cohort == target_cohort &
            ct_effects$period == target_period, ]
        if (nrow(match_row) == 0L) return(NA_real_)
        return(match_row$att[1L])
      }

    }, error = function(e) {
      NA_real_
    })
  }

  if (isTRUE(parallel) && .can_use_parallel()) {
    ri_stats <- unlist(
      .parallel_lapply(seq_len(reps), .ri_stag_one_rep),
      use.names = FALSE
    )
  } else {
    ri_stats <- vapply(seq_len(reps), .ri_stag_one_rep, numeric(1))
  }

  # -- Phase 3: Post-processing --

  ri_valid     <- ri_stats[!is.na(ri_stats)]
  n_valid      <- length(ri_valid)
  n_failed     <- reps - n_valid
  failure_rate <- n_failed / reps

  # Minimum valid threshold
  min_valid <- max(50L, as.integer(0.1 * reps))

  if (n_valid < min_valid) {
    stop_lwdid(
      sprintf(
        paste0("Staggered RI: only %d valid replications ",
               "(minimum %d required)"),
        n_valid, min_valid
      ),
      class = "lwdid_numerical"
    )
  }

  if (failure_rate > 0.1) {
    warn_lwdid(
      sprintf(
        paste0("Staggered RI: %d/%d replications failed ",
               "(%.1f%%); results may be unreliable"),
        n_failed, reps, failure_rate * 100
      ),
      class = "lwdid_numerical"
    )
  }

  # Two-sided p-value as the simple randomization tail proportion.
  ri_pvalue <- sum(abs(ri_valid) >= abs(observed_att)) / n_valid

  # Return 11-field result list
  list(
    ri_pvalue       = ri_pvalue,
    ri_distribution = ri_valid,
    obs_att         = observed_att,
    target          = target,
    target_cohort   = target_cohort,
    target_period   = target_period,
    n_valid         = n_valid,
    n_failed        = n_failed,
    failure_rate    = failure_rate,
    reps            = reps,
    method          = method
  )
}
