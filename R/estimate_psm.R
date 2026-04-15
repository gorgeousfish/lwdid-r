# ============================================================================
# estimate_psm.R — PSM ATT Estimator
# ============================================================================
# Implements Propensity Score Matching ATT with Abadie-Imbens SE.
#
# Paper references:
#   - lw2025 Procedure 3.1 Step 2 / Procedure 4.1 Step 3
#   - Abadie & Imbens (2006): matching estimator, SE, asymptotic normality
#   - Rosenbaum & Rubin (1983): propensity score theorem
# ============================================================================

# ============================================================================
# .validate_psm_inputs() — Input validation for PSM estimator
# ============================================================================
#' @keywords internal
.validate_psm_inputs <- function(data, y, d, propensity_controls,
                                  n_neighbors, match_order) {
  # 1. y column existence
  if (!y %in% names(data)) {
    stop_lwdid(sprintf("Column '%s' not found in data.", y),
               class = "lwdid_missing_column", column = y)
  }
  # 2. d column existence
  if (!d %in% names(data)) {
    stop_lwdid(sprintf("Column '%s' not found in data.", d),
               class = "lwdid_missing_column", column = d)
  }
  # 3. propensity_controls columns existence
  missing_cols <- setdiff(propensity_controls, names(data))
  if (length(missing_cols) > 0L) {
    stop_lwdid(
      sprintf("Propensity control columns not found: %s",
              paste(missing_cols, collapse = ", ")),
      class = "lwdid_missing_column", column = missing_cols)
  }
  # 4. n_neighbors >= 1
  if (!is.numeric(n_neighbors) || length(n_neighbors) != 1L ||
      n_neighbors < 1L) {
    stop_lwdid("n_neighbors must be an integer >= 1.",
               class = "lwdid_invalid_param", param = "n_neighbors")
  }
  # 5. D binary check
  d_vals <- unique(data[[d]])
  if (!all(d_vals %in% c(0, 1))) {
    stop_lwdid(
      sprintf("Treatment indicator '%s' must be binary (0/1).", d),
      class = "lwdid_invalid_data")
  }
  # 6. match_order validation
  valid_orders <- c("data", "random", "largest", "smallest")
  if (!match_order %in% valid_orders) {
    stop_lwdid(
      sprintf("Invalid match_order '%s'. Must be one of: %s",
              match_order, paste(valid_orders, collapse = ", ")),
      class = "lwdid_invalid_param", param = "match_order")
  }
  invisible(NULL)
}

# ============================================================================
# .nearest_neighbor_match() — Brute-force nearest neighbor PS matching
# ============================================================================
# Abadie & Imbens (2006) Section 2: M(i) = argmin_{j:D_j=0} |p_i - p_j|
# Receives already-reordered PS vectors; order restoration in estimate_psm().
#' @keywords internal
.nearest_neighbor_match <- function(pscores_treat, pscores_control,
                                     n_neighbors = 1L,
                                     with_replacement = TRUE,
                                     caliper = NULL) {
  n_treat  <- length(pscores_treat)
  n_control <- length(pscores_control)

  matched_control_ids <- vector("list", n_treat)
  match_counts        <- integer(n_treat)
  n_dropped           <- 0L

  # Without-replacement tracking (mimics Python set())
  used_controls <- integer(0)

  for (i in seq_len(n_treat)) {
    ps_i <- pscores_treat[i]
    distances <- abs(pscores_control - ps_i)

    # Without replacement: mark used controls as unavailable
    if (!with_replacement && length(used_controls) > 0L) {
      available <- !seq_len(n_control) %in% used_controls
      if (!any(available)) {
        # No controls left
        matched_control_ids[[i]] <- integer(0)
        match_counts[i] <- 0L
        n_dropped <- n_dropped + 1L
        next
      }
      distances[!available] <- Inf
    }

    # Caliper check: any within range?
    if (!is.null(caliper)) {
      if (!any(distances <= caliper)) {
        matched_control_ids[[i]] <- integer(0)
        match_counts[i] <- 0L
        n_dropped <- n_dropped + 1L
        next
      }
    }

    # Find k nearest neighbors
    k <- min(n_neighbors, n_control)
    nearest_indices <- order(distances)[seq_len(k)]

    # Caliper re-validation on selected neighbors
    if (!is.null(caliper)) {
      nearest_indices <- nearest_indices[distances[nearest_indices] <= caliper]
    }

    if (length(nearest_indices) == 0L) {
      matched_control_ids[[i]] <- integer(0)
      match_counts[i] <- 0L
      n_dropped <- n_dropped + 1L
      next
    }

    # Record matches
    matched_control_ids[[i]] <- nearest_indices
    match_counts[i] <- length(nearest_indices)

    # Without replacement: mark as used
    if (!with_replacement) {
      used_controls <- c(used_controls, nearest_indices)
    }
  }

  list(
    matched_control_ids = matched_control_ids,
    match_counts        = match_counts,
    n_dropped           = n_dropped
  )
}

# ============================================================================
# .compute_psm_se_abadie_imbens() — Simplified Neyman-type SE
# ============================================================================
# SE = sqrt(var(tau_i) / N_valid), where var() uses ddof=1 (R default).
# CI: att +/- qnorm(1-alpha/2) * se (normal distribution).
# NOT the full A&I (2006) Theorem 4 formula (no K_M(j), no sigma^2(X)).
#' @keywords internal
.compute_psm_se_abadie_imbens <- function(Y_treat, Y_control,
                                           matched_control_ids,
                                           att, alpha) {
  individual_effects <- numeric(0)
  for (i in seq_along(matched_control_ids)) {
    matches <- matched_control_ids[[i]]
    if (length(matches) > 0L) {
      effect_i <- Y_treat[i] - mean(Y_control[matches])
      individual_effects <- c(individual_effects, effect_i)
    }
  }

  n_valid <- length(individual_effects)
  if (n_valid < 2L) {
    return(list(se = NaN, ci_lower = NaN, ci_upper = NaN))
  }

  # var() uses ddof=1 by default, matching Python np.var(ddof=1)
  var_att <- stats::var(individual_effects) / n_valid
  se <- sqrt(var_att)

  z_crit <- stats::qnorm(1 - alpha / 2)
  ci_lower <- att - z_crit * se
  ci_upper <- att + z_crit * se

  list(se = se, ci_lower = ci_lower, ci_upper = ci_upper)
}

# ============================================================================
# .compute_psm_se_reuse_adjusted() — Reuse-adjusted normal CI
# ============================================================================
# Starts from the simplified matched-effect variance and adds a control-side
# covariance inflation term when with-replacement matching reuses controls
# across treated units.
#' @keywords internal
.compute_psm_se_reuse_adjusted <- function(Y_treat, Y_control,
                                            matched_control_ids,
                                            att, alpha) {
  individual_effects <- numeric(0)
  control_weights <- numeric(length(Y_control))

  individual_effects <- numeric(0)
  baseline_term <- 0
  treated_pos <- 0L
  for (matches in matched_control_ids) {
    treated_pos <- treated_pos + 1L
    if (length(matches) > 0L) {
      k <- length(matches)
      effect_i <- Y_treat[treated_pos] - mean(Y_control[matches])
      individual_effects <- c(individual_effects, effect_i)
    }
  }

  n_valid <- length(individual_effects)
  if (n_valid < 2L) {
    return(list(se = NaN, ci_lower = NaN, ci_upper = NaN))
  }

  treated_pos <- 0L
  for (matches in matched_control_ids) {
    treated_pos <- treated_pos + 1L
    if (length(matches) > 0L) {
      k <- length(matches)
      weight_increment <- 1 / (n_valid * k)
      control_weights[matches] <- control_weights[matches] + weight_increment
      baseline_term <- baseline_term + 1 / (n_valid^2 * k)
    }
  }

  base_var_att <- stats::var(individual_effects) / n_valid
  control_var <- stats::var(Y_control)
  if (!is.finite(control_var)) {
    control_var <- 0
  }

  reuse_term <- sum(control_weights^2)
  reuse_inflation <- control_var * max(0, reuse_term - baseline_term)
  se <- sqrt(base_var_att + reuse_inflation)

  z_crit <- stats::qnorm(1 - alpha / 2)
  ci_lower <- att - z_crit * se
  ci_upper <- att + z_crit * se

  list(se = se, ci_lower = ci_lower, ci_upper = ci_upper)
}

# ============================================================================
# .compute_psm_se_bootstrap() — Full-process resampling bootstrap
# ============================================================================
# Each bootstrap iteration: resample -> PS estimation -> matching -> ATT
# Implemented via recursive call to estimate_psm(se_method="abadie_imbens")
# CI: percentile method via quantile()
# SE: sd(att_valid) with ddof=1 (matches Python np.std(ddof=1))
#' @keywords internal
.compute_psm_se_bootstrap <- function(data, y, d, propensity_controls,
                                       n_neighbors, with_replacement,
                                       caliper, caliper_scale,
                                       trim_threshold, n_bootstrap,
                                       seed, alpha, match_order) {
  n <- nrow(data)
  att_boots <- rep(NA_real_, n_bootstrap)

  if (!is.null(seed)) set.seed(seed)

  for (b in seq_len(n_bootstrap)) {
    idx    <- sample.int(n, n, replace = TRUE)
    data_b <- data[idx, , drop = FALSE]

    # Check minimum group sizes in bootstrap sample
    D_b       <- as.integer(data_b[[d]])
    n_treat_b <- sum(D_b == 1L)
    n_ctrl_b  <- sum(D_b == 0L)
    if (n_treat_b < 1L || n_ctrl_b < n_neighbors) next

    att_boots[b] <- tryCatch({
      res_b <- estimate_psm(
        data = data_b, y = y, d = d,
        propensity_controls = propensity_controls,
        n_neighbors = n_neighbors,
        with_replacement = with_replacement,
        caliper = caliper, caliper_scale = caliper_scale,
        trim_threshold = trim_threshold,
        se_method = "abadie_imbens",
        seed = NULL, alpha = alpha,
        match_order = match_order
      )
      res_b$att
    }, error = function(e) NA_real_)
  }

  att_valid    <- att_boots[!is.na(att_boots)]
  n_valid_boot <- length(att_valid)
  success_rate <- n_valid_boot / n_bootstrap

  # Convergence warning if success rate < 50%
  if (success_rate < 0.5) {
    warn_lwdid(
      sprintf("Bootstrap success rate %.1f%% (< 50%%).", success_rate * 100),
      class = "lwdid_convergence")
  }

  # Error if fewer than 10 valid bootstrap samples
  if (n_valid_boot < 10L) {
    stop_lwdid(
      sprintf("Only %d valid bootstrap samples (< 10). Cannot estimate SE.",
              n_valid_boot),
      class = "lwdid_estimation_failed")
  }

  se       <- stats::sd(att_valid)
  ci_lower <- as.numeric(stats::quantile(att_valid, alpha / 2))
  ci_upper <- as.numeric(stats::quantile(att_valid, 1 - alpha / 2))

  list(se = se, ci_lower = ci_lower, ci_upper = ci_upper)
}

# ============================================================================
# estimate_psm() — PSM ATT Main Function
# ============================================================================
# lw2025 Procedure 3.1 Step 2 / Procedure 4.1 Step 3 (matching option)
# Abadie & Imbens (2006) Equations 1-2: counterfactual construction & ATT
#' Estimate PSM-ATT (Propensity Score Matching)
#'
#' @param data data.frame, cross-sectional data (one row per unit)
#' @param y character, outcome variable name
#' @param d character, treatment indicator name (1=treated, 0=control)
#' @param propensity_controls character vector, PS model covariates
#' @param n_neighbors integer, number of nearest neighbors (default 1)
#' @param with_replacement logical, match with replacement (default TRUE)
#' @param caliper numeric or NULL, caliper width
#' @param caliper_scale character, "sd" or "absolute"
#' @param trim_threshold numeric, PS trimming threshold (default 0.01)
#' @param se_method character, "abadie_imbens", "reuse_adjusted", or "bootstrap"
#' @param n_bootstrap integer, bootstrap replications (default 200)
#' @param seed integer or NULL, random seed
#' @param alpha numeric, significance level (default 0.05)
#' @param match_order character, matching order
#'
#' @return list with 16 fields
#' @export
estimate_psm <- function(data, y, d, propensity_controls,
                         n_neighbors = 1L, with_replacement = TRUE,
                         caliper = NULL, caliper_scale = "sd",
                         trim_threshold = 0.01, se_method = "abadie_imbens",
                         n_bootstrap = 200L, seed = NULL,
                         alpha = 0.05, match_order = "data") {

  # ---- Step 1: Input validation ----
  .validate_psm_inputs(data, y, d, propensity_controls,
                       n_neighbors, match_order)

  # Validate se_method
  if (!se_method %in% c("abadie_imbens", "reuse_adjusted", "bootstrap")) {
    stop_lwdid(
      sprintf(
        "Invalid se_method '%s'. Use 'abadie_imbens', 'reuse_adjusted', or 'bootstrap'.",
              se_method),
      class = "lwdid_invalid_param", param = "se_method")
  }

  # ---- Step 2: Missing value exclusion ----
  all_vars <- unique(c(y, d, propensity_controls))
  complete_mask <- stats::complete.cases(data[, all_vars, drop = FALSE])
  data_clean <- data[complete_mask, , drop = FALSE]

  # ---- Step 3: Sample size checks ----
  D <- as.integer(data_clean[[d]])
  Y <- as.numeric(data_clean[[y]])
  n_treated <- sum(D == 1L)
  n_control <- sum(D == 0L)

  if (n_treated < 1L) {
    stop_lwdid("No treated units found. Cannot perform PSM.",
               class = "lwdid_no_treated")
  }
  if (n_control < n_neighbors) {
    stop_lwdid(
      sprintf("Control sample size (%d) < n_neighbors (%d).",
              n_control, n_neighbors),
      class = "lwdid_insufficient_data")
  }

  # Small sample warning
  if (n_treated < 5L) {
    warn_lwdid(
      sprintf("Small treated sample (n=%d). Estimates may be unstable.",
              n_treated),
      class = "lwdid_small_sample")
  }

  # ---- Step 4: Propensity score estimation ----
  ps_result <- estimate_propensity_score(
    data_clean, d, propensity_controls, trim_threshold)
  ps <- ps_result$propensity_scores

  # ---- Step 5: Extreme PS diagnostic ----
  prop_low  <- mean(ps < 0.05)
  prop_high <- mean(ps > 0.95)
  if (prop_low > 0.1 || prop_high > 0.1) {
    warn_lwdid(
      sprintf("Extreme PS: %.1f%% < 0.05, %.1f%% > 0.95. Overlap concern.",
              prop_low * 100, prop_high * 100),
      class = "lwdid_overlap")
  }

  # ---- Step 6: Match order ----
  treat_idx <- which(D == 1L)
  ctrl_idx  <- which(D == 0L)
  ps_treat  <- ps[treat_idx]
  ps_ctrl   <- ps[ctrl_idx]
  Y_treat   <- Y[treat_idx]
  Y_ctrl    <- Y[ctrl_idx]

  if (!is.null(seed)) set.seed(seed)

  order_idx <- switch(match_order,
    "data"     = seq_len(n_treated),
    "random"   = sample.int(n_treated),
    "largest"  = order(ps_treat, decreasing = TRUE),
    "smallest" = order(ps_treat, decreasing = FALSE)
  )

  ps_treat_ordered <- ps_treat[order_idx]

  # ---- Step 7: Caliper conversion ----
  actual_caliper <- caliper
  if (!is.null(caliper)) {
    if (caliper_scale == "sd") {
      # ddof=0 population std, matching Python np.std()
      ps_mean <- mean(ps)
      ps_sd_pop <- sqrt(sum((ps - ps_mean)^2) / length(ps))
      actual_caliper <- caliper * ps_sd_pop
    }
    # caliper_scale == "absolute": no conversion needed
  }

  # ---- Step 8: Matching ----
  match_result <- .nearest_neighbor_match(
    ps_treat_ordered, ps_ctrl,
    n_neighbors = n_neighbors,
    with_replacement = with_replacement,
    caliper = actual_caliper
  )

  # ---- Step 9: Restore original order ----
  # match_result is in order_idx order; map back to original treat order
  matched_ids_original    <- vector("list", n_treated)
  match_counts_original   <- integer(n_treated)
  for (j in seq_len(n_treated)) {
    orig_pos <- order_idx[j]
    matched_ids_original[[orig_pos]]  <- match_result$matched_control_ids[[j]]
    match_counts_original[orig_pos]   <- match_result$match_counts[j]
  }
  n_dropped <- match_result$n_dropped

  # ---- Step 10: ATT calculation ----
  # A&I (2006) Eq 1-2: tau_i = Y_i - mean(Y_control[M(i)])
  att_individual <- numeric(0)
  for (i in seq_len(n_treated)) {
    matches <- matched_ids_original[[i]]
    if (length(matches) > 0L) {
      effect_i <- Y_treat[i] - mean(Y_ctrl[matches])
      att_individual <- c(att_individual, effect_i)
    }
  }

  n_valid <- length(att_individual)

  # ---- Step 11: Matching diagnostics ----
  if (n_valid == 0L) {
    stop_lwdid("All treated units failed to match. Cannot estimate ATT.",
               class = "lwdid_estimation_failed")
  }

  match_success_rate <- n_valid / n_treated
  if (match_success_rate < 0.5) {
    warn_lwdid(
      sprintf("Low match success rate (%.1f%%). Results may be unreliable.",
              match_success_rate * 100),
      class = "lwdid_overlap")
  }

  att <- mean(att_individual)

  # n_matched: unique matched control units
  all_matched <- unlist(matched_ids_original)
  n_matched <- length(unique(all_matched[!is.na(all_matched)]))

  # ---- Step 12: SE calculation ----
  if (se_method == "abadie_imbens") {
    se_result <- .compute_psm_se_abadie_imbens(
      Y_treat, Y_ctrl, matched_ids_original, att, alpha)
  } else if (se_method == "reuse_adjusted") {
    se_result <- .compute_psm_se_reuse_adjusted(
      Y_treat, Y_ctrl, matched_ids_original, att, alpha)
  } else {
    se_result <- .compute_psm_se_bootstrap(
      data_clean, y, d, propensity_controls,
      n_neighbors, with_replacement,
      caliper, caliper_scale,
      trim_threshold, n_bootstrap,
      seed, alpha, match_order)
  }

  se       <- se_result$se
  ci_lower <- se_result$ci_lower
  ci_upper <- se_result$ci_upper

  # ---- Step 13: Normal inference ----
  if (is.finite(se) && se > 0) {
    t_stat <- att / se
    pvalue <- 2 * (1 - stats::pnorm(abs(t_stat)))
  } else {
    t_stat <- NaN
    pvalue <- NaN
  }

  # ---- Step 14: Return 16-field list ----
  list(
    att                = att,
    se                 = se,
    ci_lower           = ci_lower,
    ci_upper           = ci_upper,
    t_stat             = t_stat,
    pvalue             = pvalue,
    propensity_scores  = ps,
    match_counts       = match_counts_original,
    matched_control_ids = matched_ids_original,
    n_treated          = n_treated,
    n_control          = n_control,
    n_matched          = n_matched,
    caliper            = actual_caliper,
    n_dropped          = n_dropped,
    match_success_rate = match_success_rate,
    estimator          = "psm"
  )
}
