# lwdid-r/R/lwdid.R

#' Lee-Wooldridge Difference-in-Differences Estimation
#'
#' @description
#' Estimate treatment effects using the Lee-Wooldridge (2025)
#' Difference-in-Differences framework. Supports both Common Timing
#' (all units treated at the same time) and Staggered Adoption
#' (units treated at different times) designs.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Capture the original call via \code{match.call()}
#'   \item Initialize and reset the global WarningRegistry
#'   \item Process backward-compatible \code{riseed} argument
#'   \item Validate all inputs via \code{\link{validate_inputs}}
#'   \item Dispatch to Common Timing or Staggered estimation
#'   \item Flush warnings and extract diagnostics
#'   \item Populate metadata (call, version)
#'   \item Optionally produce a plot
#' }
#'
#' @param data data.frame or data.table. Panel data in long format.
#'   Each row represents one unit-period observation.
#' @param y character(1). Name of the outcome variable column.
#' @param ivar character(1). Name of the unit identifier column.
#'   Required -- no default value.
#' @param tvar character(1) or character(2). Name(s) of the time variable
#'   column(s). Length 1 for annual data (e.g., \code{"year"});
#'   length 2 for quarterly data as \code{c(year_var, quarter_var)}.
#'   Required -- no default value.
#' @param d character(1) or NULL. Name of the treatment indicator column
#'   (Common Timing mode). Must be time-invariant within each unit.
#' @param post character(1) or NULL. Name of the post-treatment indicator
#'   column (Common Timing mode).
#' @param gvar character(1) or NULL. Name of the cohort variable column
#'   (Staggered mode). Indicates the first treatment period for each unit.
#'   Never-treated units encoded as NA, 0, or Inf.
#' @param control_group character(1). Control group strategy for Staggered
#'   mode. One of \code{"not_yet_treated"}, \code{"never_treated"},
#'   \code{"all_others"}, \code{"auto"}. Default \code{"not_yet_treated"}.
#' @param aggregate character(1). Aggregation type for Staggered mode.
#'   One of \code{"none"}, \code{"cohort"}, \code{"overall"}.
#'   Default \code{"cohort"}.
#' @param event_time_range Numeric vector of length 2 or NULL. Used when
#'   \code{aggregate = "event_time"} to keep only event times in
#'   \code{c(min_e, max_e)}.
#' @param df_strategy Character(1). Degrees-of-freedom aggregation strategy for
#'   event-time summaries. One of \code{"conservative"}, \code{"weighted"}, or
#'   \code{"fallback"}. Default \code{"conservative"}.
#' @param rolling character(1). Transformation method. One of
#'   \code{"demean"}, \code{"detrend"}, \code{"demeanq"}, \code{"detrendq"}.
#'   Case-insensitive. Default \code{"demean"}.

#' @param estimator character(1). Estimator type. One of \code{"ra"},
#'   \code{"ipw"}, \code{"ipwra"}, \code{"psm"}. Default \code{"ra"}.
#' @param controls character vector or NULL. Names of control variable columns.
#' @param ps_controls character vector or NULL. Names of propensity score
#'   model control variable columns.
#' @param vce character(1) or NULL. Variance-covariance estimator type.
#'   One of \code{NULL}, \code{"hc0"}-\code{"hc4"}, \code{"robust"},
#'   \code{"cluster"}, \code{"bootstrap"}.
#' @param cluster_var character(1) or NULL. Name of the cluster variable
#'   column. Required when \code{vce} is \code{"cluster"} or
#'   \code{"bootstrap"}.
#' @param alpha numeric(1). Significance level for confidence intervals.
#'   Default 0.05.
#' @param ri logical(1). Whether to perform randomization inference.
#'   Default FALSE.
#' @param rireps integer(1). Number of RI replications. Default 1000L.
#' @param seed integer(1) or NULL. Random seed for reproducibility.
#' @param ri_method character(1). RI method. One of \code{"bootstrap"},
#'   \code{"permutation"}. Default \code{"bootstrap"}.
#' @param season_var character(1) or NULL. Name of the seasonal indicator
#'   column. Required for \code{demeanq}/\code{detrendq} when \code{tvar}
#'   is length 1.
#' @param quarter character(1) or NULL. Backward-compatible alias for
#'   \code{season_var}.
#' @param Q integer(1). Number of seasons per cycle. Default 4L.
#' @param auto_detect_frequency logical(1). Whether to auto-detect data
#'   frequency. Default FALSE.
#' @param include_pretreatment logical(1). Whether to include pre-treatment
#'   period estimates. Default FALSE.
#' @param pretreatment_test logical(1). Whether to perform parallel trends
#'   test. Default TRUE.
#' @param pretreatment_alpha numeric(1). Significance level for
#'   pre-treatment test. Default 0.05.
#' @param exclude_pre_periods integer(1). Number of pre-treatment periods
#'   to exclude from estimation. Default 0L.
#' @param wcb_reps integer(1). Number of wild cluster bootstrap replications.
#'   Default 999L.
#' @param wcb_type character(1). Wild bootstrap weight type. One of
#'   \code{"rademacher"}, \code{"mammen"}, or \code{"webb"}.
#' @param wcb_seed integer(1) or NULL. Optional seed passed to
#'   \code{\link{wild_cluster_bootstrap}}.
#' @param wcb_restricted_model character(1). Restricted model passed to
#'   \code{\link{wild_cluster_bootstrap}}. One of
#'   \code{"with_controls"} or \code{"intercept_only"}.
#' @param auto_wcb logical(1). When \code{TRUE}, \code{vce="cluster"} with
#'   fewer than 20 clusters automatically upgrades to wild cluster bootstrap.
#' @param use_fwildclusterboot logical(1). Whether to request the
#'   \code{fwildclusterboot} adapter before falling back to the native R path.
#' @param trim_threshold numeric(1). Propensity score trimming threshold.
#'   Default 0.01.
#' @param trim_method character(1). Trimming method. One of \code{"clip"},
#'   \code{"drop"}. Default \code{"clip"}.
#' @param n_neighbors integer(1). Number of PSM neighbors. Default 1L.
#' @param caliper numeric(1) or NULL. PSM caliper distance.
#' @param caliper_scale character(1). Caliper scale. One of \code{"sd"},
#'   \code{"absolute"}. Default \code{"sd"}.
#' @param with_replacement logical(1). Whether PSM uses replacement.
#'   Default TRUE.
#' @param match_order character(1). PSM match order. One of \code{"data"},
#'   \code{"random"}, \code{"largest"}, \code{"smallest"}.
#'   Default \code{"data"}.
#' @param se_method character(1) or NULL. Standard error method.
#'   For IPW/IPWRA: \code{"analytical"} or \code{"bootstrap"}.
#'   For PSM: \code{"abadie_imbens"} or \code{"bootstrap"}.
#'   NULL uses estimator default. Ignored for RA.
#' @param boot_reps integer(1). Number of bootstrap replications.
#'   Default 200L.
#' @param return_diagnostics logical(1). Whether to return propensity score
#'   diagnostics. Default FALSE.
#' @param parallel logical(1). Whether to enable parallel computation for
#'   staggered (g,r) estimation loop. Requires \pkg{future} and
#'   \pkg{future.apply}. Default FALSE.
#' @param n_cores integer(1) or NULL. Number of parallel workers. NULL
#'   auto-detects via \code{future::availableCores() - 1}. Ignored when
#'   \code{parallel = FALSE}.
#' @param balanced_panel character(1). Panel balance handling. One of
#'   \code{"warn"}, \code{"error"}, \code{"ignore"}. Default \code{"warn"}.
#' @param graph logical(1). Whether to produce a plot after estimation.
#'   Default FALSE.
#' @param gid NULL or scalar. Unit ID to highlight in the plot.
#' @param graph_options list or NULL. Additional options passed to the
#'   plot method.
#' @param verbose character(1). Verbosity level. One of \code{"quiet"},
#'   \code{"default"}, \code{"verbose"}. Default \code{"default"}.
#' @param att_gt logical(1). R extension parameter. Default FALSE.
#' @param ... Additional arguments. Currently supports \code{riseed}
#'   (deprecated alias for \code{seed}). Unknown arguments trigger a
#'   warning.
#'
#' @return An S3 object of class \code{\link{lwdid_result}} containing
#'   estimation results, diagnostics, and metadata.
#'
#' @details
#' When \code{estimator = "ra"} (default), inference uses the exact
#' t-distribution (Lee & Wooldridge 2026, Equation 2.10). For
#' \code{estimator = "ipw"}, \code{"ipwra"}, or \code{"psm"}, inference
#' uses the asymptotic normal distribution based on:
#' \itemize{
#'   \item IPW: Lunceford & Davidian (2004) asymptotic normality
#'   \item IPWRA: Cattaneo (2010) asymptotic normality
#'   \item PSM: Abadie & Imbens (2006) asymptotic normality
#' }
#' This distinction affects confidence interval construction: RA uses
#' \code{qt(1 - alpha/2, df)} critical values, while IPW/IPWRA/PSM
#' use \code{qnorm(1 - alpha/2)}.
#'
#' @examples
#' \dontrun{
#' # Common Timing design
#' result <- lwdid(data = panel_df, y = "outcome", ivar = "unit_id",
#'                 tvar = "year", d = "treated", post = "post_period")
#' print(result)
#' summary(result)
#'
#' # Staggered Adoption design
#' result <- lwdid(data = panel_df, y = "outcome", ivar = "unit_id",
#'                 tvar = "year", gvar = "first_treat_year")
#' print(result)
#' }
#'
#' @export
#' @seealso
#' \code{\link{validate_inputs}} for input validation,
#' \code{\link{new_lwdid_result}} for result object structure,
#' \code{\link{print.lwdid_result}} for printing results,
#' \code{\link{plot.lwdid_result}} for plotting,
#' \code{\link{plot_event_study.lwdid_result}} for event study plots

lwdid <- function(
  data,
  y,
  ivar,
  tvar,
  d = NULL,
  post = NULL,
  gvar = NULL,
  control_group = "not_yet_treated",
  aggregate = "cohort",
  event_time_range = NULL,
  df_strategy = "conservative",
  rolling = "demean",
  estimator = "ra",
  controls = NULL,
  ps_controls = NULL,
  vce = NULL,
  cluster_var = NULL,
  alpha = 0.05,
  ri = FALSE,
  rireps = 1000L,
  seed = NULL,
  ri_method = c("bootstrap", "permutation"),
  season_var = NULL,
  quarter = NULL,
  Q = 4L,
  auto_detect_frequency = FALSE,
  include_pretreatment = FALSE,
  pretreatment_test = TRUE,
  pretreatment_alpha = 0.05,
  exclude_pre_periods = 0L,
  wcb_reps = 999L,
  wcb_type = c("rademacher", "mammen", "webb"),
  wcb_seed = NULL,
  wcb_restricted_model = "with_controls",
  auto_wcb = TRUE,
  use_fwildclusterboot = TRUE,
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
  parallel = FALSE,
  n_cores = NULL,
  balanced_panel = "warn",
  graph = FALSE,
  gid = NULL,
  graph_options = NULL,
  verbose = "default",
  att_gt = FALSE,
  ...
) {

  # ---------------------------------------------------------------------------
  # Step 1: Capture the original call (Task E1-06.2)
  # ---------------------------------------------------------------------------
  cl <- match.call()

  # ---------------------------------------------------------------------------
  # Step 2: Initialize and reset the global WarningRegistry (Task E1-06.3)
  # ---------------------------------------------------------------------------
  registry <- .lwdid_env$warning_registry
  registry$clear()

  # ---------------------------------------------------------------------------
  # Step 3: Process backward-compatible kwargs (Task E1-06.5)
  # ---------------------------------------------------------------------------
  dots <- list(...)
  known_kwargs <- "riseed"

  # 3a: riseed backward compatibility
  if (!is.null(dots$riseed)) {
    if (is.null(seed)) {
      riseed_val <- dots$riseed
      if (is.numeric(riseed_val) || is.integer(riseed_val)) {
        seed <- as.integer(riseed_val)
      } else if (is.character(riseed_val)) {
        seed_parsed <- suppressWarnings(as.integer(riseed_val))
        if (is.na(seed_parsed)) {
          stop_lwdid(
            message = paste0(
              "riseed must be an integer or integer-convertible string, got '",
              riseed_val, "'. Use seed parameter for explicit control."
            ),
            class = "lwdid_invalid_parameter",
            param = "riseed",
            value = riseed_val
          )
        }
        seed <- seed_parsed
      } else {
        stop_lwdid(
          message = paste0(
            "riseed must be an integer, got ", class(riseed_val)[1],
            ". Use seed parameter for explicit control."
          ),
          class = "lwdid_invalid_parameter",
          param = "riseed",
          value = riseed_val
        )
      }
    }
    dots$riseed <- NULL
  }

  # 3b: Unknown kwargs warning
  unknown_kwargs <- setdiff(names(dots), known_kwargs)
  if (length(unknown_kwargs) > 0) {
    registry$register(
      category = "lwdid_data",
      message = sprintf(
        "Unknown argument(s) ignored: %s. Valid extra arguments: %s.",
        paste(unknown_kwargs, collapse = ", "),
        paste(known_kwargs, collapse = ", "))
    )
  }

  # ---------------------------------------------------------------------------
  # Step 3c: Validate new RI / pretreatment parameters (Task E7-06.1)
  # ---------------------------------------------------------------------------
  ri_method <- match.arg(ri_method)
  wcb_type <- match.arg(wcb_type)

  stopifnot(is.logical(ri), length(ri) == 1L)
  stopifnot(is.logical(include_pretreatment), length(include_pretreatment) == 1L)
  stopifnot(is.logical(pretreatment_test), length(pretreatment_test) == 1L)
  stopifnot(is.numeric(rireps), length(rireps) == 1L, rireps >= 1L)
  stopifnot(is.null(seed) || (is.numeric(seed) && length(seed) == 1L && seed >= 1L))
  stopifnot(is.numeric(pretreatment_alpha), length(pretreatment_alpha) == 1L,
            pretreatment_alpha > 0, pretreatment_alpha < 1)

  if (!is.numeric(exclude_pre_periods) || length(exclude_pre_periods) != 1L ||
      exclude_pre_periods < 0 || exclude_pre_periods != as.integer(exclude_pre_periods)) {
    stop_lwdid("exclude_pre_periods must be a non-negative integer",
               class = "lwdid_invalid_param")
  }
  if (!is.numeric(wcb_reps) || length(wcb_reps) != 1L ||
      wcb_reps < 1 || wcb_reps != as.integer(wcb_reps)) {
    stop_lwdid(
      "wcb_reps must be a positive integer",
      class = "lwdid_invalid_parameter",
      param = "wcb_reps",
      value = wcb_reps,
      allowed = "positive integer(1)"
    )
  }
  if (!is.null(wcb_seed) &&
      (!is.numeric(wcb_seed) || length(wcb_seed) != 1L ||
       wcb_seed < 1 || wcb_seed != as.integer(wcb_seed))) {
    stop_lwdid(
      "wcb_seed must be NULL or a positive integer",
      class = "lwdid_invalid_parameter",
      param = "wcb_seed",
      value = wcb_seed,
      allowed = "NULL or positive integer(1)"
    )
  }
  if (!wcb_restricted_model %in% c("with_controls", "intercept_only")) {
    stop_lwdid(
      "wcb_restricted_model must be 'with_controls' or 'intercept_only'",
      class = "lwdid_invalid_parameter",
      param = "wcb_restricted_model",
      value = wcb_restricted_model,
      allowed = c("with_controls", "intercept_only")
    )
  }
  if (!is.logical(auto_wcb) || length(auto_wcb) != 1L || is.na(auto_wcb)) {
    stop_lwdid(
      "auto_wcb must be TRUE or FALSE",
      class = "lwdid_invalid_parameter",
      param = "auto_wcb",
      value = auto_wcb,
      allowed = c(TRUE, FALSE)
    )
  }
  if (!is.logical(use_fwildclusterboot) || length(use_fwildclusterboot) != 1L ||
      is.na(use_fwildclusterboot)) {
    stop_lwdid(
      "use_fwildclusterboot must be TRUE or FALSE",
      class = "lwdid_invalid_parameter",
      param = "use_fwildclusterboot",
      value = use_fwildclusterboot,
      allowed = c(TRUE, FALSE)
    )
  }

  rireps <- as.integer(rireps)
  exclude_pre_periods <- as.integer(exclude_pre_periods)
  wcb_reps <- as.integer(wcb_reps)

  season_var_effective <- if (!is.null(season_var) || !is.null(quarter)) {
    .resolve_season_var(season_var = season_var, quarter = quarter)
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # Step 4: Validate all inputs (Task E1-06.4)
  # ---------------------------------------------------------------------------
  validated <- validate_inputs(
    data = data, y = y, ivar = ivar, tvar = tvar,
    d = d, post = post, gvar = gvar,
    control_group = control_group, aggregate = aggregate,
    rolling = rolling, estimator = estimator,
    controls = controls, ps_controls = ps_controls,
    vce = vce, cluster_var = cluster_var, alpha = alpha,
    ri = ri, rireps = rireps, seed = seed, ri_method = ri_method,
    season_var = season_var_effective, Q = Q,
    auto_detect_frequency = auto_detect_frequency,
    include_pretreatment = include_pretreatment,
    pretreatment_test = pretreatment_test,
    pretreatment_alpha = pretreatment_alpha,
    exclude_pre_periods = exclude_pre_periods,
    trim_threshold = trim_threshold, trim_method = trim_method,
    n_neighbors = n_neighbors,
    caliper = caliper, caliper_scale = caliper_scale,
    with_replacement = with_replacement,
    match_order = match_order, se_method = se_method,
    boot_reps = boot_reps,
    return_diagnostics = return_diagnostics,
    balanced_panel = balanced_panel,
    graph = graph, gid = gid, graph_options = graph_options,
    verbose = verbose, att_gt = att_gt,
    registry = registry
  )

  # Add E5-05 params to validated_params (not yet in validate_inputs)
  validated$validated_params$event_time_range <- event_time_range
  validated$validated_params$df_strategy <- df_strategy
  validated$validated_params$quarter <- quarter
  validated$validated_params$wcb_reps <- wcb_reps
  validated$validated_params$wcb_type <- wcb_type
  validated$validated_params$wcb_seed <- wcb_seed
  validated$validated_params$wcb_restricted_model <- wcb_restricted_model
  validated$validated_params$auto_wcb <- auto_wcb
  validated$validated_params$use_fwildclusterboot <- use_fwildclusterboot
  validated$validated_params$parallel <- parallel
  validated$validated_params$n_cores <- n_cores

  # ---------------------------------------------------------------------------
  # Step 5: Mode dispatch (Task E1-06.6)
  # ---------------------------------------------------------------------------
  result <- if (validated$mode == "staggered") {
    .estimate_staggered(validated, registry = registry)
  } else {
    .estimate_common_timing(validated, registry = registry)
  }

  # ---------------------------------------------------------------------------
  # Step 5a: Pre-treatment Integration (Task E7-06.3)
  # ---------------------------------------------------------------------------
  if (isTRUE(include_pretreatment)) {
    vp <- validated$validated_params
    pre_effects <- NULL

    # Convert tpost1 from tindex scale to original time scale.
    # validated$tpost1 is in tindex (1-based sequential), but
    # estimate_pre_treatment_common() compares against data[[tvar]]
    # which is in original scale (e.g., 2005). Formula:
    #   original_tpost1 = T_min + tpost1_tindex - 1
    tpost1_orig <- validated$T_min + (validated$tpost1 %||% result$tpost1) - 1L

    tryCatch({
      if (!isTRUE(result$is_staggered)) {
        # CT pretreatment (Task E7-06.3.1)
        pre_effects <- estimate_pre_treatment_common(
          data = data,
          y = y, ivar = ivar, d = vp$d %||% d,
          tvar = tvar,
          tpost1 = tpost1_orig,
          rolling = rolling,
          controls = controls,
          estimator = vp$estimator %||% "ra",
          vce = vce,
          cluster_var = cluster_var,
          alpha = alpha,
          exclude_pre_periods = exclude_pre_periods
        )
      } else {
        # Staggered pretreatment (Task E7-06.3.2)
        pre_effects <- estimate_pre_treatment_staggered(
          data = data,
          y = y, ivar = ivar, tvar = tvar,
          gvar = gvar,
          rolling = rolling,
          estimator = vp$estimator %||% "ra",
          controls = controls,
          control_group = vp$control_group %||% "not_yet_treated",
          vce = vce,
          cluster_var = cluster_var,
          alpha = alpha
        )
      }
    }, error = function(e) {
      warn_lwdid(
        sprintf("Pre-treatment estimation failed: %s",
                conditionMessage(e)),
        class = "lwdid_data"
      )
    })

    result$att_pre_treatment <- pre_effects

    # Joint test (Task E7-06.3.4)
    if (isTRUE(pretreatment_test) && !is.null(pre_effects) &&
        is.data.frame(pre_effects) && nrow(pre_effects) > 0L) {
      tryCatch({
        result$parallel_trends_test <- test_parallel_trends_joint(
          pre_treatment_effects = pre_effects,
          alpha = pretreatment_alpha,
          test_type = "f"
        )
      }, error = function(e) {
        warn_lwdid(
          sprintf("Parallel trends joint test failed: %s",
                  conditionMessage(e)),
          class = "lwdid_data"
        )
      })
    }

  }

  # Always store exclude_pre_periods in result (even when include_pretreatment=FALSE)
  result$exclude_pre_periods <- as.integer(exclude_pre_periods)

  # ---------------------------------------------------------------------------
  # Step 5b: Wild cluster bootstrap integration (Epic 9)
  # ---------------------------------------------------------------------------
  if (!isTRUE(validated$mode == "staggered")) {
    result$wcb_auto_triggered <- FALSE

    vce_lower <- if (is.null(validated$validated_params$vce)) {
      NULL
    } else {
      tolower(validated$validated_params$vce)
    }

    use_wcb <- FALSE
    auto_wcb_triggered <- FALSE
    n_clusters_wcb <- NA_integer_

    if (identical(vce_lower, "bootstrap")) {
      use_wcb <- TRUE
    } else if (identical(vce_lower, "cluster") && isTRUE(auto_wcb)) {
      if (!is.null(result$.wcb_data) &&
          !is.null(cluster_var) &&
          cluster_var %in% names(result$.wcb_data)) {
        n_clusters_wcb <- length(unique(result$.wcb_data[[cluster_var]]))
      } else if (!is.null(result$n_clusters) && !is.na(result$n_clusters)) {
        n_clusters_wcb <- as.integer(result$n_clusters)
      }

      if (is.finite(n_clusters_wcb) && n_clusters_wcb < 20L) {
        use_wcb <- TRUE
        auto_wcb_triggered <- TRUE
        message(sprintf(
          "Number of clusters G=%d < 20, auto-enabling Wild Cluster Bootstrap",
          as.integer(n_clusters_wcb)
        ))
      }
    }

    if (isTRUE(use_wcb)) {
      wcb_data <- result$.wcb_data
      wcb_cluster_var <- result$.wcb_cluster_var %||% cluster_var
      if (is.null(wcb_data) || is.null(wcb_cluster_var)) {
        stop_lwdid(
          "Wild cluster bootstrap integration requires prepared common-timing cross-section data and cluster identifiers.",
          class = "lwdid_internal_error"
        )
      }

      use_fwildclusterboot_actual <- isTRUE(use_fwildclusterboot)
      if (use_fwildclusterboot_actual && !.fwildclusterboot_available()) {
        message(
          paste0(
            "fwildclusterboot not installed, using pure R implementation. Install via: ",
            "install.packages('fwildclusterboot', repos='https://s3alfisc.r-universe.dev')"
          )
        )
        use_fwildclusterboot_actual <- FALSE
      }

      wcb_result <- wild_cluster_bootstrap(
        data = wcb_data,
        y_transformed = result$.wcb_y_transformed %||% "y_trans_summary",
        d = result$.wcb_d %||% "d",
        cluster_var = wcb_cluster_var,
        controls = controls,
        n_bootstrap = wcb_reps,
        weight_type = wcb_type,
        alpha = alpha,
        seed = wcb_seed,
        impose_null = TRUE,
        restricted_model = wcb_restricted_model,
        use_fwildclusterboot = use_fwildclusterboot_actual
      )

      result$pvalue_wcb <- wcb_result$pvalue
      result$ci_lower_wcb <- wcb_result$ci_lower
      result$ci_upper_wcb <- wcb_result$ci_upper
      result$se_wcb <- wcb_result$se_bootstrap
      result$wcb_details <- wcb_result
      result$analytic_inference <- list(
        se_att = result$se_att,
        pvalue = result$pvalue,
        ci_lower = result$ci_lower,
        ci_upper = result$ci_upper,
        t_stat = result$t_stat,
        df_inference = result$df_inference,
        vce_type = result$vce_type
      )
      result$se_att <- result$se_wcb
      result$pvalue <- result$pvalue_wcb
      result$ci_lower <- result$ci_lower_wcb
      result$ci_upper <- result$ci_upper_wcb
      result$t_stat <- wcb_result$t_stat_original %||% NA_real_
      result$vce_type <- "wild_cluster_bootstrap"
      result$wcb_auto_triggered <- auto_wcb_triggered
    }
  }

  # ---------------------------------------------------------------------------
  # Step 5c: Randomization Inference integration (Task E7-06.2)
  # ---------------------------------------------------------------------------
  is_staggered <- (validated$mode == "staggered")

  if (isTRUE(ri)) {
    # --- E7-06.2.1: Seed management ---
    actual_seed <- if (!is.null(seed)) {
      as.integer(seed)
    } else {
      sample.int(.Machine$integer.max, 1L)
    }

    # Initialize RI result fields
    ri_pvalue <- NULL
    ri_distribution <- NULL
    ri_valid <- NULL
    ri_failed <- NULL
    ri_error <- NULL
    ri_target <- NULL

    if (!is_staggered) {
      # --- E7-06.2.2: Common Timing RI ---
      # Extract firstpost cross-section transformed data from result
      y_trans_ri <- result$.ri_y_trans
      d_ri <- result$.ri_d
      x_ri <- result$.ri_x

      tryCatch({
        ri_result <- randomization_inference(
          y_trans = y_trans_ri,
          d = d_ri,
          x = x_ri,
          reps = rireps,
          seed = actual_seed,
          method = ri_method
        )
        ri_pvalue <- ri_result$ri_pvalue
        ri_distribution <- ri_result$ri_distribution
        ri_valid <- ri_result$n_valid
        ri_failed <- ri_result$n_failed
      }, error = function(e) {
        # --- E7-06.2.3: CT RI failure degradation ---
        warn_lwdid(
          sprintf(
            "Randomization inference failed: %s. ATT estimation results are still valid.",
            conditionMessage(e)
          ),
          class = "lwdid_numerical"
        )
        ri_error <<- conditionMessage(e)
        ri_valid <<- 0L
        ri_failed <<- -1L
      })

    } else {
      # --- E7-06.2.4-E7-06.2.8: Staggered RI ---
      att_overall_val <- result$att_overall
      aggregate_val <- validated$validated_params$aggregate
      ri_observed <- NULL
      target_cohort_ri <- NULL
      target_period_ri <- NULL

      # --- E7-06.2.4: Nested if-else target mapping ---
      if (aggregate_val == "overall") {
        # --- E7-06.2.5: Overall branch ---
        if (!is.null(att_overall_val) && !is.na(att_overall_val)) {
          ri_target <- "overall"
          ri_observed <- att_overall_val
          target_cohort_ri <- NULL
          target_period_ri <- NULL
        } else {
          ri_error <- "Overall ATT is NULL or NA"
          ri_target <- "overall"
          ri_valid <- 0L
          ri_failed <- -1L
        }
      } else if (aggregate_val == "cohort") {
        # --- E7-06.2.6: Cohort branch ---
        cohort_effs <- result$cohort_effects
        if (!is.null(cohort_effs) && length(cohort_effs) > 0L) {
          # Filter for valid (non-NA) ATT values
          valid_ce <- Filter(function(e) {
            !is.null(e$att) && is.finite(e$att)
          }, cohort_effs)
          if (length(valid_ce) > 0L) {
            ri_target <- "cohort"
            ri_observed <- valid_ce[[1L]]$att
            target_cohort_ri <- as.integer(valid_ce[[1L]]$cohort)
            target_period_ri <- NULL
          } else {
            ri_error <- "No valid cohort ATT estimates available"
            ri_target <- "cohort"
            ri_valid <- 0L
            ri_failed <- -1L
          }
        } else {
          ri_error <- "No cohort ATT estimates available"
          ri_target <- "cohort"
          ri_valid <- 0L
          ri_failed <- -1L
        }
      } else {
        # --- E7-06.2.7: Cohort-time branch ---
        ri_target <- "cohort_time"
        ct_effects <- result$att_by_cohort_time
        if (!is.null(ct_effects) && nrow(ct_effects) > 0L) {
          valid_rows <- ct_effects[!is.na(ct_effects$att), , drop = FALSE]
          if (nrow(valid_rows) > 0L) {
            ri_observed <- valid_rows$att[1L]
            target_cohort_ri <- as.integer(valid_rows$cohort[1L])
            target_period_ri <- as.integer(valid_rows$period[1L])
          } else {
            ri_error <- "No valid cohort-time ATT estimates available"
            ri_valid <- 0L
            ri_failed <- -1L
          }
        } else {
          ri_error <- "No available effect estimates"
          ri_valid <- 0L
          ri_failed <- -1L
        }
      }

      # --- E7-06.2.8: Execute Staggered RI (only if ri_observed determined) ---
      if (!is.null(ri_observed)) {
        tryCatch({
          ri_result <- ri_staggered(
            data = validated$data,
            y = validated$validated_params$depvar,
            ivar = validated$validated_params$ivar,
            tvar = validated$validated_params$tvar,
            gvar = validated$validated_params$gvar,
            observed_att = ri_observed,
            target = ri_target,
            target_cohort = target_cohort_ri,
            target_period = target_period_ri,
            rolling = validated$validated_params$rolling,
            controls = validated$validated_params$controls,
            vce = validated$validated_params$vce,
            cluster_var = validated$validated_params$cluster_var,
            reps = rireps,
            seed = actual_seed,
            method = ri_method
          )
          ri_pvalue <- ri_result$ri_pvalue
          ri_distribution <- ri_result$ri_distribution
          ri_valid <- ri_result$n_valid
          ri_failed <- ri_result$n_failed
        }, error = function(e) {
          warn_lwdid(
            sprintf(
              "Randomization inference failed: %s",
              conditionMessage(e)
            ),
            class = "lwdid_numerical"
          )
          ri_error <<- conditionMessage(e)
          ri_valid <<- 0L
          ri_failed <<- -1L
        })
      }
    }

    # --- E7-06.2.9: Store RI results into result object ---
    result$ri_pvalue <- ri_pvalue
    result$ri_distribution <- ri_distribution
    result$ri_seed <- actual_seed
    result$ri_method <- ri_method
    result$rireps <- rireps
    result$ri_valid <- ri_valid
    result$ri_failed <- ri_failed
    result$ri_error <- ri_error
    result$ri_target <- ri_target
  }

  # Clean up internal RI data fields (not part of public API)
  result$.ri_y_trans <- NULL
  result$.ri_d <- NULL
  result$.ri_x <- NULL

  if (isTRUE(validated$validated_params$return_diagnostics)) {
    result <- .attach_mainline_diagnostics(result, validated)
  }

  # ---------------------------------------------------------------------------
  # Step 6: Registry flush and diagnostics extraction (Task E1-06.7)
  # ---------------------------------------------------------------------------
  registry$flush(verbose = validated$validated_params$verbose)
  result$warning_diagnostics <- registry$get_diagnostics()
  result$warnings_log <- registry$get_log()

  # ---------------------------------------------------------------------------
  # Step 7: Metadata population (Task E1-06.8)
  # ---------------------------------------------------------------------------
  result$call <- cl
  result$lwdid_version <- as.character(utils::packageVersion("lwdid"))

  # ---------------------------------------------------------------------------
  # Step 8: Graph processing (Task E1-06.9)
  # ---------------------------------------------------------------------------
  if (isTRUE(validated$validated_params$graph)) {
    tryCatch(
      plot(result, gid = validated$validated_params$gid,
           graph_options = validated$validated_params$graph_options),
      error = function(e) {
        registry$register(
          category = "lwdid_data",
          message = sprintf(
            "Plotting failed: %s: %s. The estimation results are unaffected.",
            class(e)[1], conditionMessage(e))
        )
      }
    )
  }

  # ---------------------------------------------------------------------------
  # Step 9: Return result
  # ---------------------------------------------------------------------------
  result
}


# =============================================================================
# Helper Functions
# =============================================================================

#' Extract unit-level control variable matrix
#'
#' Extract time-invariant control variables for each unit. If controls
#' vary over time within the pre-treatment window, compute pre-period
#' means and issue a warning.
#'
#' @param dt data.table, panel data.
#' @param ivar character, unit identifier column name.
#' @param tvar character, time variable column name.
#' @param controls character vector, control variable column names.
#' @param S integer, treatment time.
#' @param exclude_pre_periods integer, number of pre-periods to exclude.
#' @return data.table with ivar column + K control variable columns
#'   (one row per unit).
#' @keywords internal
.extract_controls <- function(dt, ivar, tvar, controls,
                              S, exclude_pre_periods) {
  pre_end <- S - 1L - exclude_pre_periods

  # Pre-period data for control extraction
  pre_dt <- dt[get(tvar) <= pre_end]

  # Check time-invariance for each control variable
  time_varying <- character(0)
  for (ctrl in controls) {
    n_unique <- pre_dt[, .(nu = data.table::uniqueN(
      stats::na.omit(get(ctrl))
    )), by = ivar]
    if (any(n_unique$nu > 1L)) {
      time_varying <- c(time_varying, ctrl)
    }
  }

  if (length(time_varying) > 0L) {
    warn_lwdid(
      sprintf(
        paste0("Control variable(s) [%s] vary over time; ",
               "pre-period means used as time-invariant proxy"),
        paste(time_varying, collapse = ", ")
      ),
      class = "lwdid_data",
      detail = "controls_time_varying",
      action_taken = "pre-period means used"
    )
  }

  # Compute per-unit control values (mean over pre-period)
  x_dt <- pre_dt[, lapply(.SD, function(v) mean(v, na.rm = TRUE)),
                 by = ivar, .SDcols = controls]

  # Ensure all units present (units with no pre-period rows get NA)
  all_units <- data.table::data.table(
    tmp_ivar = unique(dt[[ivar]])
  )
  data.table::setnames(all_units, "tmp_ivar", ivar)
  x_dt <- merge(all_units, x_dt, by = ivar, all.x = TRUE)

  x_dt
}

# =============================================================================
# Estimation Functions
# =============================================================================

#' Complete Common Timing estimation pipeline
#'
#' Implements the full Common Timing estimation flow: data
#' preparation, rolling transformation, post-period mean aggregation
#' (lw2026 equation 2.13), NA filtering, RA estimation,
#' period-specific effects, and result construction.
#'
#' Cluster pre-validation is performed before transformation:
#' column existence, NA check, and nesting within treatment groups
#' are verified on the original panel data (see lw2026 Section 8.2
#' for clustering discussion).
#'
#' @param validated List returned by \code{\link{validate_inputs}}.
#' @param registry WarningRegistry instance for deferred warning
#'   emission.
#'
#' @references Lee, S. and Wooldridge, J.M. (2025, 2026).
#'   Section 8.2 (clustering).
#'
#' @return An S3 object of class \code{\link{lwdid_result}}.
#' @keywords internal
.estimate_common_timing <- function(validated, registry = NULL) {

  # =========================================================================
  # Step 1a: Parameter extraction and defensive copy
  # =========================================================================
  params <- validated$validated_params
  dt <- data.table::copy(validated$data)
  S <- validated$tpost1

  y <- params$depvar
  ivar <- params$ivar
  tvar <- params$tvar
  d_var <- params$d
  rolling <- params$rolling
  alpha <- params$alpha
  exclude_pre_periods <- params$exclude_pre_periods
  controls <- params$controls
  estimator_controls <- unique(c(controls, params$ps_controls))
  if (length(estimator_controls) == 0L) {
    estimator_controls <- NULL
  }
  vce <- params$vce
  vce_for_ols <- if (!is.null(vce) && identical(tolower(vce), "bootstrap")) {
    "cluster"
  } else {
    vce
  }
  cluster_var <- params$cluster_var

  # =========================================================================
  # Step 1b: Treatment indicator and control variable extraction
  # =========================================================================
  # Unit-level treatment indicator (time-invariant, one value per unit)
  unit_d <- dt[, .(d = data.table::first(get(d_var))), by = ivar]

  # Control variable matrix (if any)
  # Pass "tindex" since S is in tindex scale
  x_dt <- NULL
  if (!is.null(estimator_controls) && length(estimator_controls) > 0L) {
    x_dt <- .extract_controls(
      dt, ivar, "tindex", estimator_controls, S, exclude_pre_periods
    )
  }

  # =========================================================================
  # Step 1c: Cluster variable pre-validation (before transformation)
  # Must run on original panel data -- after transformation each unit
  # has only one row, so nesting is trivially satisfied.
  # =========================================================================
  cluster_vals <- NULL
  cluster_count <- NA_integer_
  allow_degenerate_wcb <- identical(tolower(vce %||% ""), "bootstrap") ||
    (identical(tolower(vce %||% ""), "cluster") && isTRUE(params$auto_wcb))
  if (!is.null(vce_for_ols) &&
      identical(tolower(vce_for_ols), "cluster") &&
      !is.null(cluster_var)) {
    # Validation 1: column existence
    if (!cluster_var %in% names(dt)) {
      stop_lwdid(
        sprintf(
          "Cluster variable '%s' not found in data. Available columns: %s",
          cluster_var,
          paste(head(names(dt), 10), collapse = ", ")
        ),
        class = "lwdid_invalid_parameter",
        param = "cluster_var", value = cluster_var
      )
    }

    # Validation 2: NA check
    if (anyNA(dt[[cluster_var]])) {
      n_na <- sum(is.na(dt[[cluster_var]]))
      stop_lwdid(
        sprintf(
          "Cluster variable '%s' contains %d NA value(s). Remove or impute before using cluster VCE.",
          cluster_var, n_na
        ),
        class = "lwdid_invalid_parameter",
        param = "cluster_var", value = sprintf("contains %d NA", n_na)
      )
    }

    # Validation 3: Nesting -- each unit must belong to exactly one cluster
    unit_cluster_counts <- tapply(
      dt[[cluster_var]], dt[[ivar]],
      function(x) length(unique(x))
    )
    violating <- names(unit_cluster_counts)[unit_cluster_counts > 1L]
    if (length(violating) > 0L) {
      examples <- head(violating, 5)
      example_details <- vapply(examples, function(uid) {
        clusters <- unique(dt[[cluster_var]][dt[[ivar]] == uid])
        sprintf("%s -> {%s}", uid, paste(clusters, collapse = ", "))
      }, character(1))
      stop_lwdid(
        sprintf(
          paste0("Cluster nesting violated: %d unit(s) belong to multiple clusters. ",
                 "Examples: %s"),
          length(violating),
          paste(example_details, collapse = "; ")
        ),
        class = "lwdid_invalid_parameter",
        param = "cluster_var",
        value = sprintf("%d units with multiple clusters", length(violating))
      )
    }
  }

  # =========================================================================
  # Step 2: Rolling transformation
  # =========================================================================
  # S (tpost1) is in tindex scale, so use "tindex" as the time variable
  # for all operations that compare against S. tindex is 1-based sequential
  # (created by .create_time_index: year - min(year) + 1).
  dt <- transform_common(dt, y, ivar, "tindex", g = S,
                         rolling = rolling,
                         exclude_pre_periods = exclude_pre_periods,
                         season_var = params$season_var,
                         Q = params$Q)

  # =========================================================================
  # Step 3a: Post-period mean aggregation -- lw2026 equation 2.13
  # =========================================================================
  summary_dt <- dt[tindex >= S,
                   .(y_trans_summary = mean(y_trans, na.rm = TRUE)),
                   by = ivar]

  # Merge treatment indicator (aligned by ivar)
  summary_dt <- merge(summary_dt, unit_d, by = ivar)

  # Merge control variables (if any, aligned by ivar)
  if (!is.null(x_dt)) {
    summary_dt <- merge(summary_dt, x_dt, by = ivar)
  }

  # =========================================================================
  # Step 3b: NA filtering
  # qr() does not handle NA -- NA propagates to all coefficients.
  # y_trans_summary is NA/NaN from two sources:
  #   (1) Units with insufficient pre-periods (y_trans all NA)
  #   (2) Units with no valid post-period Y observations
  # R: is.na(NaN) returns TRUE, handling both sources uniformly.
  # =========================================================================
  valid <- !is.na(summary_dt$y_trans_summary)
  n_excluded <- sum(!valid)
  if (n_excluded > 0L) {
    warn_lwdid(
      sprintf(
        "%d unit(s) excluded due to NA transformed outcomes",
        n_excluded
      ),
      class = "lwdid_data",
      detail = "units_excluded_na_ytrans"
    )
  }
  summary_valid <- summary_dt[valid]
  if (nrow(summary_valid) < 3L) {
    stop_lwdid(
      sprintf(
        paste0("Fewer than 3 valid units after NA filtering ",
               "(%d remaining)"),
        nrow(summary_valid)
      ),
      class = "lwdid_insufficient_data"
    )
  }

  # =========================================================================
  # Step 3c: Extract regression inputs
  # Alignment guaranteed by merge(by = ivar) above.
  # =========================================================================
  y_summary <- summary_valid$y_trans_summary
  d_summary <- summary_valid$d
  x_summary <- if (!is.null(controls) && length(controls) > 0L) {
    as.matrix(summary_valid[, controls, with = FALSE])
  } else {
    NULL
  }
  est_covariates <- if (!is.null(estimator_controls) &&
                        length(estimator_controls) > 0L) {
    as.data.frame(summary_valid[, estimator_controls, with = FALSE])
  } else {
    NULL
  }

  # Extract cluster values at unit level (one per unit, post-merge)
  cluster_vals <- NULL
  if (!is.null(vce_for_ols) &&
      identical(tolower(vce_for_ols), "cluster") &&
      !is.null(cluster_var)) {
    # After merge, summary_valid has one row per unit.
    # If cluster_var already equals ivar (or has already been carried
    # forward), avoid re-merging the same column and creating duplicates.
    if (!cluster_var %in% names(summary_valid)) {
      unit_cluster <- dt[, .(
        .cluster_value = data.table::first(get(cluster_var))
      ), by = ivar]
      data.table::setnames(unit_cluster, ".cluster_value", cluster_var)
      summary_valid <- merge(summary_valid, unit_cluster, by = ivar, all.x = TRUE)
    }
    cluster_vals <- summary_valid[[cluster_var]]
    cluster_count <- length(unique(cluster_vals))
  }

  if (isTRUE(allow_degenerate_wcb) &&
      is.finite(cluster_count) &&
      cluster_count < 2L) {
    vce_for_ols <- NULL
  }

  # =========================================================================
  # Step 3d: Call estimator via dispatch_estimator()
  # For estimator="ra", dispatch_estimator routes to estimate_ra_common
  # which uses vector interface. For IPW/IPWRA/PSM, constructs a
  # data.frame with named columns for the data.frame-based estimators.
  # =========================================================================
  vp <- validated$validated_params
  estimator <- vp$estimator %||% "ra"

  if (estimator == "ra") {
    # RA path: use vector interface directly (backward compatible)
    att_result <- estimate_ra_common(
      y_summary, d_summary, x_summary,
      vce = vce_for_ols, cluster = cluster_vals, alpha = alpha
    )
    att_result$estimator <- "ra"
    att_result$inference_dist <- "t"
  } else {
    # IPW/IPWRA/PSM path: construct data.frame for dispatch_estimator
    est_df <- data.frame(
      .y_outcome = y_summary,
      .d_treat = d_summary
    )
    if (!is.null(est_covariates)) {
      est_df <- cbind(est_df, est_covariates)
    }
    att_result <- dispatch_estimator(
      data = est_df,
      y = ".y_outcome",
      d = ".d_treat",
      controls = controls,
      ps_controls = vp$ps_controls,
      estimator = estimator,
      vce = vce_for_ols,
      cluster_var = NULL,
      alpha = alpha,
      trim_threshold = vp$trim_threshold %||% 0.01,
      trim_method = vp$trim_method %||% "clip",
      n_neighbors = vp$n_neighbors %||% 1L,
      caliper = vp$caliper,
      caliper_scale = vp$caliper_scale %||% "sd",
      with_replacement = vp$with_replacement %||% TRUE,
      match_order = vp$match_order %||% "data",
      se_method = vp$se_method,
      boot_reps = vp$boot_reps %||% 200L,
      seed = vp$seed,
      return_diagnostics = vp$return_diagnostics %||% FALSE
    )
  }

  # =========================================================================
  # Step 4: Period-specific effects
  # Calls estimate_period_effects() (Story E2-05).
  # If not yet implemented, use placeholder.
  # =========================================================================
  periods <- sort(unique(dt[tindex >= S, tindex]))

  att_by_period <- tryCatch(
    estimate_period_effects(
      dt_transformed = dt,
      y_trans_col = "y_trans",
      d_col = d_var,
      tvar = "tindex",
      x = controls,
      periods = periods,
      vce = vce_for_ols,
      cluster_var = if (!is.null(vce_for_ols) && tolower(vce_for_ols) == "cluster")
                      cluster_var else NULL,
      alpha = alpha,
      include_pretreatment = FALSE
    ),
    error = function(e) {
      # If estimate_period_effects() not yet implemented, return NULL.
      # Match both English and localized (e.g. Chinese) error messages.
      msg <- conditionMessage(e)
      if (grepl("could not find function", msg, fixed = TRUE) ||
          grepl("estimate_period_effects", msg, fixed = TRUE)) {
        NULL
      } else {
        stop(e)
      }
    }
  )

  # Convert tindex-based period column back to original time scale
  # tindex = original_year - T_min + 1, so original_year = tindex + T_min - 1
  if (!is.null(att_by_period) &&
      is.data.frame(att_by_period) &&
      nrow(att_by_period) > 0L) {
    T_min_val <- validated$T_min
    if (!is.null(T_min_val) && is.finite(T_min_val)) {
      att_by_period$period <- att_by_period$tindex + T_min_val - 1L
    }
  }

  if (!is.null(att_by_period) &&
      is.data.frame(att_by_period) &&
      nrow(att_by_period) > 0L &&
      length(tvar) == 2L &&
      rolling %in% c("demeanq", "detrendq")) {
    period_labels <- .build_ct_public_period_labels(dt, tvar)
    att_by_period <- .decorate_ct_public_att_by_period(
      att_by_period = att_by_period,
      att_result = att_result,
      period_labels = period_labels
    )
  }

  # =========================================================================
  # Step 5: (VCE now integrated via estimate_ra_common)
  # =========================================================================

  # =========================================================================
  # Step 6: Construct lwdid_result object
  # =========================================================================
  result <- new_lwdid_result(
    # --- Core estimates ---
    att = att_result$att,
    se_att = att_result$se,
    t_stat = att_result$t_stat,
    pvalue = att_result$pvalue,
    ci_lower = att_result$ci_lower,
    ci_upper = att_result$ci_upper,
    df_resid = att_result$df,
    df_inference = att_result$df,
    # --- Sample information ---
    nobs = att_result$n,
    n_treated = att_result$n_treated,
    n_control = att_result$n_control,
    # Public result K is the last pre-treatment period index, not the
    # estimator-side control-count field returned by estimate_ra_common().
    K = validated$K %||% att_result$K,
    tpost1 = S,
    # --- Model specification ---
    depvar = y,
    rolling = rolling,
    estimator = params$estimator,
    method = "common_timing",
    vce_type = att_result$vce_type,
    cluster_var = params$cluster_var,
    n_clusters = att_result$n_clusters,
    alpha = alpha,
    is_staggered = FALSE,
    controls_used = !is.null(controls) && length(controls) > 0L,
    controls = controls %||% character(0),
    include_pretreatment = params$include_pretreatment,
    control_group = NA_character_,
    control_group_used = NA_character_,
    # --- Full regression results ---
    params = att_result$params,
    bse = att_result$se,
    vcov_matrix = att_result$vcov,
    resid = att_result$resid,
    data = NULL,
    # --- Period-specific ---
    att_by_period = att_by_period,
    att_pre_treatment = NULL,
    parallel_trends_test = NULL,
    # --- Staggered-specific (all NULL for Common Timing) ---
    aggregate = NULL,
    cohorts = NULL, cohort_sizes = NULL,
    n_never_treated = NULL,
    att_by_cohort = NULL, att_by_cohort_time = NULL,
    cohort_effects = NULL,
    att_overall = NULL, se_overall = NULL,
    ci_overall_lower = NULL, ci_overall_upper = NULL,
    t_stat_overall = NULL, pvalue_overall = NULL,
    event_time_effects = NULL, cohort_weights = NULL,
    # --- RI (NULL for now, populated by lwdid() main function) ---
    ri_pvalue = NULL, ri_distribution = NULL,
    ri_seed = NULL, rireps = NULL, ri_method = NULL,
    ri_valid = NULL, ri_failed = NULL,
    ri_error = NULL, ri_target = NULL,
    # --- Diagnostics ---
    diagnostics = list(
      controls_tier = att_result$controls_tier
    ),
    warning_diagnostics = list(),
    propensity_scores = NULL,
    matched_data = NULL, n_matched = NULL,
    match_rate = NULL,
    weights_cv = NULL, warnings_log = list(),
    # --- Metadata (populated by lwdid() main function) ---
    call = NULL, lwdid_version = NULL,
    ivar = ivar,
    tvar = tvar,
    is_quarterly = validated$is_quarterly %||% FALSE
  )

  # Store RI data for CT RI integration (Task E7-06.2)
  result$.ri_data <- list(
    y_trans = y_summary,
    d = d_summary,
    x = x_summary
  )

  # Store cross-section data for CT RI (used by lwdid() main function)
  result$.ri_y_trans <- y_summary
  result$.ri_d <- d_summary
  result$.ri_x <- x_summary
  result$.wcb_data <- as.data.frame(summary_valid)
  result$.wcb_y_transformed <- "y_trans_summary"
  result$.wcb_d <- "d"
  result$.wcb_cluster_var <- cluster_var

  result
}

.format_ct_public_period_component <- function(value) {
  if (length(value) != 1L || is.na(value)) {
    return(NA_character_)
  }

  if (is.numeric(value)) {
    if (value == as.integer(value)) {
      return(as.character(as.integer(value)))
    }
    return(as.character(value))
  }

  as.character(value)
}

.build_ct_public_period_labels <- function(dt, tvar) {
  year_var <- tvar[[1L]]
  season_var <- tvar[[2L]]

  label_dt <- dt[
    ,
    .(
      period = {
        year_val <- data.table::first(get(year_var))
        season_val <- data.table::first(get(season_var))
        year_str <- .format_ct_public_period_component(year_val)
        season_str <- .format_ct_public_period_component(season_val)
        if (is.na(year_str) || is.na(season_str)) {
          paste0("T", data.table::first(tindex))
        } else {
          paste0(year_str, "q", season_str)
        }
      }
    ),
    by = tindex
  ]

  setNames(label_dt$period, as.character(label_dt$tindex))
}

.decorate_ct_public_att_by_period <- function(att_by_period,
                                              att_result,
                                              period_labels) {
  period_rows <- att_by_period
  period_rows$tindex <- as.character(period_rows$tindex)
  period_rows$period <- vapply(
    period_rows$tindex,
    function(idx) {
      period_labels[[idx]] %||% idx
    },
    character(1)
  )

  average_row <- data.frame(
    tindex = "-",
    period = "average",
    att = unname(as.numeric(att_result$att)),
    se = unname(as.numeric(att_result$se)),
    t_stat = unname(as.numeric(att_result$t_stat)),
    pvalue = unname(as.numeric(att_result$pvalue)),
    ci_lower = unname(as.numeric(att_result$ci_lower)),
    ci_upper = unname(as.numeric(att_result$ci_upper)),
    n_obs = as.integer(att_result$n),
    n_treated = as.integer(att_result$n_treated),
    n_control = as.integer(att_result$n_control),
    df = if (is.null(att_result$df)) NA_integer_ else as.integer(att_result$df),
    vce_type = att_result$vce_type %||% NA_character_,
    n_clusters = if (is.null(att_result$n_clusters)) NA_integer_ else as.integer(att_result$n_clusters),
    controls_tier = att_result$controls_tier %||% NA_character_,
    stringsAsFactors = FALSE
  )

  out <- rbind(average_row, period_rows)
  rownames(out) <- NULL
  out
}


#' Aggregate (g,r) effects to cohort-level for non-RA estimators
#'
#' Computes cohort ATTs as equal-weight averages of period-specific
#' ATT(g,r) estimates. Used when estimator is IPW/IPWRA/PSM where
#' the time-averaged RA regression (equations 7.9-7.10) is not
#' applicable.
#'
#' @param cohort_time_effects data.frame with columns: cohort, period,
#'   att, se, n_treated, n_control, df_inference
#' @param cohorts integer vector of cohort values
#' @param alpha numeric significance level
#' @return list of cohort effect results matching aggregate_to_cohort format
#' @keywords internal
.aggregate_gr_to_cohort <- function(cohort_time_effects, cohorts, alpha = 0.05) {
  results <- list()
  for (g in cohorts) {
    gr_g <- cohort_time_effects[
      cohort_time_effects$cohort == g & !is.na(cohort_time_effects$att), ,
      drop = FALSE
    ]
    if (nrow(gr_g) == 0L) next
    n_periods <- nrow(gr_g)
    att_g <- mean(gr_g$att)
    # Conservative SE: assume independence across r within cohort g
    se_vals <- gr_g$se[is.finite(gr_g$se)]
    se_g <- if (length(se_vals) > 0L) sqrt(sum(se_vals^2)) / n_periods else NA_real_
    # Use median df from (g,r) effects
    df_vals <- gr_g$df_inference[is.finite(gr_g$df_inference)]
    df_g <- if (length(df_vals) > 0L) {
      as.integer(stats::median(df_vals))
    } else {
      NA_integer_
    }
    # Compute inference
    t_stat_g <- if (is.finite(se_g) && se_g > 0) att_g / se_g else NA_real_
    pvalue_g <- if (!is.na(df_g) && df_g > 0L && !is.na(t_stat_g)) {
      2 * stats::pt(abs(t_stat_g), df = df_g, lower.tail = FALSE)
    } else {
      NA_real_
    }
    t_crit <- if (!is.na(df_g) && df_g > 0L) {
      stats::qt(1 - alpha / 2, df = df_g)
    } else {
      stats::qnorm(1 - alpha / 2)
    }
    ci_lower_g <- att_g - t_crit * se_g
    ci_upper_g <- att_g + t_crit * se_g
    n_treat <- if ("n_treated" %in% names(gr_g)) {
      max(gr_g$n_treated, na.rm = TRUE)
    } else { NA_integer_ }
    n_ctrl <- if ("n_control" %in% names(gr_g)) {
      max(gr_g$n_control, na.rm = TRUE)
    } else { NA_integer_ }
    results[[length(results) + 1L]] <- list(
      cohort = as.integer(g),
      att = att_g,
      se = se_g,
      ci_lower = ci_lower_g,
      ci_upper = ci_upper_g,
      t_stat = t_stat_g,
      pvalue = pvalue_g,
      n_periods = n_periods,
      n_units = n_treat,
      n_control = n_ctrl,
      df_resid = df_g,
      df_inference = df_g,
      K = NA_integer_
    )
  }
  results
}

.estimate_staggered <- function(validated, registry = NULL) {
  # --- Extract convenience aliases from validated ---
  dt <- validated$data
  y <- validated$validated_params$depvar
  ivar <- validated$validated_params$ivar
  tvar <- validated$validated_params$tvar
  gvar <- validated$validated_params$gvar
  rolling <- validated$validated_params$rolling
  controls <- validated$validated_params$controls
  vce <- validated$validated_params$vce
  cluster_var <- validated$validated_params$cluster_var
  alpha <- validated$validated_params$alpha
  aggregate <- validated$validated_params$aggregate
  exclude_pre_periods <- as.integer(
    validated$validated_params$exclude_pre_periods
  )
  vp <- validated$validated_params
  estimator <- vp$estimator %||% "ra"
  resolved_cg <- validated$control_group_resolved
  cohorts <- validated$cohorts
  n_nt <- validated$n_nt
  T_max <- validated$T_max
  event_time_range <- vp$event_time_range
  df_strategy <- vp$df_strategy %||% "conservative"
  verbose <- vp$verbose %||% "default"
  verbose_flag <- !identical(verbose, "silent")
  parallel <- vp$parallel %||% FALSE
  n_cores <- vp$n_cores

  # --- Step 0: Parameter validation (E5-05.1) ---
  valid_aggregates <- c("none", "cohort", "overall", "event_time")
  if (!aggregate %in% valid_aggregates) {
    stop_lwdid(
      sprintf("Invalid aggregate value '%s'. Must be one of: %s",
              aggregate, paste(valid_aggregates, collapse = ", ")),
      class = "lwdid_invalid_input"
    )
  }

  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop_lwdid(
      sprintf("alpha must be in (0, 1), got: %s", deparse(alpha)),
      class = "lwdid_invalid_input"
    )
  }

  if (!is.null(event_time_range)) {
    if (!is.numeric(event_time_range) || length(event_time_range) != 2L ||
        any(is.na(event_time_range)) ||
        event_time_range[1] > event_time_range[2]) {
      stop_lwdid(
        sprintf("event_time_range must be length-2 numeric with lower <= upper, got: %s",
                deparse(event_time_range)),
        class = "lwdid_invalid_input"
      )
    }
  }

  valid_df_strategies <- c("conservative", "weighted", "fallback")
  if (!df_strategy %in% valid_df_strategies) {
    stop_lwdid(
      sprintf("Invalid df_strategy '%s'. Must be one of: %s",
              df_strategy, paste(valid_df_strategies, collapse = ", ")),
      class = "lwdid_invalid_input"
    )
  }

  # --- Small NT sample warning (R enhancement: broader scope) ---
  if (resolved_cg == "never_treated" && n_nt < 2L) {
    warn_lwdid(
      sprintf(
        paste("Never-treated unit count is very small",
              "(N=%d); inference results may be unreliable.",
              "Recommend N_NT >= 2"),
        n_nt),
      class = "lwdid_data",
      detail = "small_nt_sample"
    )
  }

  # --- Step 3: Precompute transform statistics ---
  pre_stats <- precompute_transforms(
    dt = dt, y = y, ivar = ivar, tvar = tvar,
    cohorts = cohorts, rolling = rolling,
    exclude_pre_periods = exclude_pre_periods
  )
  valid_cohorts <- sort(as.integer(names(pre_stats)))

  if (length(valid_cohorts) == 0L) {
    stop_lwdid(
      paste("All cohorts have insufficient pre-treatment",
            "periods; cannot perform transformation"),
      class = "lwdid_insufficient_data"
    )
  }
  if (length(valid_cohorts) < length(cohorts)) {
    n_dropped <- length(cohorts) - length(valid_cohorts)
    message(sprintf(
      paste("[lwdid] %d cohort(s) skipped due to insufficient",
            "pre-periods; %d valid cohort(s) remaining"),
      n_dropped, length(valid_cohorts)
    ))
  }

  # --- Step 4: (g,r) effect estimation ---
  cohort_time_effects <- estimate_staggered_effects(
    dt = dt, y = y, ivar = ivar, tvar = tvar,
    gvar = gvar, rolling = rolling,
    control_group = resolved_cg, controls = controls,
    vce = vce, cluster_var = cluster_var,
    alpha = alpha, pre_stats = pre_stats,
    estimator = vp$estimator %||% "ra",
    ps_controls = vp$ps_controls,
    trim_threshold = vp$trim_threshold %||% 0.01,
    trim_method = vp$trim_method %||% "clip",
    n_neighbors = vp$n_neighbors %||% 1L,
    caliper = vp$caliper,
    caliper_scale = vp$caliper_scale %||% "sd",
    with_replacement = vp$with_replacement %||% TRUE,
    match_order = vp$match_order %||% "data",
    se_method = vp$se_method,
    boot_reps = vp$boot_reps %||% 200L,
    seed = vp$seed,
    return_diagnostics = vp$return_diagnostics %||% FALSE,
    parallel = parallel,
    n_cores = n_cores
  )

  # --- Step 1: Control group constraint check (FATAL-004) ---
  control_group_used <- resolved_cg
  if (aggregate %in% c("cohort", "overall")) {
    # Need NT units for cohort/overall aggregation
    if (n_nt == 0L) {
      stop_lwdid(
        sprintf(
          paste0("No never-treated units found. aggregate='%s' requires ",
                 "never-treated control group. Consider using ",
                 "aggregate='none' or aggregate='event_time'."),
          aggregate),
        class = "lwdid_no_never_treated"
      )
    }
    if (n_nt < 2L) {
      warn_lwdid(
        sprintf(
          paste0("Only %d never-treated unit(s) found for aggregate='%s'. ",
                 "Inference may be unreliable. Recommend N_NT >= 2."),
          n_nt, aggregate),
        class = "lwdid_small_sample",
        detail = "small_nt_for_aggregation"
      )
    }
    # Auto-switch control group to never_treated
    if (resolved_cg != "never_treated") {
      warn_lwdid(
        sprintf(
          paste0("Control group auto-switched from '%s' to 'never_treated' ",
                 "for aggregate='%s'. Cohort/overall aggregation requires ",
                 "never-treated units as the control group."),
          resolved_cg, aggregate),
        class = "lwdid_control_group_switch"
      )
      control_group_used <- "never_treated"
      # Re-run (g,r) estimation with never_treated control group
      cohort_time_effects <- estimate_staggered_effects(
        dt = dt, y = y, ivar = ivar, tvar = tvar,
        gvar = gvar, rolling = rolling,
        control_group = "never_treated", controls = controls,
        vce = vce, cluster_var = cluster_var,
        alpha = alpha, pre_stats = pre_stats,
        estimator = vp$estimator %||% "ra",
        ps_controls = vp$ps_controls,
        trim_threshold = vp$trim_threshold %||% 0.01,
        trim_method = vp$trim_method %||% "clip",
        n_neighbors = vp$n_neighbors %||% 1L,
        caliper = vp$caliper,
        caliper_scale = vp$caliper_scale %||% "sd",
        with_replacement = vp$with_replacement %||% TRUE,
        match_order = vp$match_order %||% "data",
        se_method = vp$se_method,
        boot_reps = vp$boot_reps %||% 200L,
        seed = vp$seed,
        return_diagnostics = vp$return_diagnostics %||% FALSE,
        parallel = parallel,
        n_cores = n_cores
      )
    }
  }
  # event_time does NOT trigger control group switch

  # --- Step 4: Cohort sizes computation ---
  cohort_sizes <- validated$cohort_sizes

  # --- Step 5: Aggregation dispatch ---
  cohort_effects <- NULL
  overall_effect <- NULL
  event_time_effects <- NULL
  att_cohort_agg <- NULL
  se_cohort_agg <- NULL
  t_stat_cohort_agg <- NULL
  pvalue_cohort_agg <- NULL
  ci_cohort_agg <- NULL

  if (aggregate %in% c("cohort", "overall")) {
    if (estimator == "ra") {
      # Step 5a: Cohort aggregation via time-averaged regression
      # (lw2026 equations 7.9-7.10) — exact for RA estimator
      cohort_effects <- aggregate_to_cohort(
        dt = dt, y = y, ivar = ivar, tvar = tvar, gvar = gvar,
        cohorts = valid_cohorts, T_max = T_max,
        pre_stats = pre_stats, rolling = rolling,
        vce = vce, cluster_var = cluster_var,
        alpha = alpha, controls = controls
      )
    } else {
      # Step 5a-alt: For non-RA estimators (IPW/IPWRA/PSM), compute
      # cohort effects as equal-weight averages of (g,r) effects.
      # aggregate_to_cohort uses RA regression which is not appropriate
      # for IPW/IPWRA/PSM where the ATT is estimated via different
      # weighting/matching procedures.
      cohort_effects <- .aggregate_gr_to_cohort(
        cohort_time_effects, valid_cohorts, alpha
      )
    }
  }

  if (aggregate == "overall") {
    # Step 5b: Overall aggregation
    overall_effect <- aggregate_to_overall(
      dt = dt, y = y, ivar = ivar, tvar = tvar, gvar = gvar,
      cohorts = valid_cohorts, T_max = T_max,
      pre_stats = pre_stats, rolling = rolling,
      vce = vce, cluster_var = cluster_var,
      alpha = alpha, controls = controls
      # TODO (Epic 9): Map transform_type for demeanq/detrendq
      # Python maps demeanq -> 'detrend' in aggregate_to_overall
      # Note: Python mapping inconsistency between cohort and overall
    )
  }

  if (aggregate == "event_time") {
    # Step 5c: Event-time aggregation
    event_time_effects <- aggregate_to_event_time(
      cohort_time_effects = cohort_time_effects,
      cohort_sizes = cohort_sizes,
      event_time_range = event_time_range,
      df_strategy = df_strategy,
      alpha = alpha,
      verbose = verbose_flag
    )
  }

  # --- Step 5a-pre: df_fallback computation ---
  n_controls_k <- if (!is.null(controls)) length(controls) else 0L
  n_total_units <- as.integer(validated$n_treated) + as.integer(n_nt)
  df_fallback <- max(1L, n_total_units - 2L - 2L * n_controls_k)

  # --- Step 5b: Cohort aggregation statistics ---
  if (aggregate == "cohort" && !is.null(cohort_effects)) {
    # Filter valid cohort effects
    valid_ce <- Filter(function(e) {
      is.finite(e$att) && is.finite(e$se)
    }, cohort_effects)

    if (length(valid_ce) > 0L) {
      ce_n_units <- vapply(valid_ce, function(e) {
        as.numeric(e$n_units)
      }, numeric(1))
      ce_att <- vapply(valid_ce, `[[`, numeric(1), "att")
      ce_se <- vapply(valid_ce, `[[`, numeric(1), "se")

      total_units <- sum(ce_n_units)
      w_g <- ce_n_units / total_units

      att_cohort_agg <- sum(w_g * ce_att)
      se_cohort_agg <- sqrt(sum(w_g^2 * ce_se^2))

      if (is.finite(se_cohort_agg) && se_cohort_agg > 0) {
        t_stat_cohort_agg <- att_cohort_agg / se_cohort_agg

        # df for cohort agg: median of cohort df_inference
        ce_df_inf <- vapply(valid_ce, function(e) {
          as.numeric(e$df_inference)
        }, numeric(1))
        valid_ce_df <- ce_df_inf[is.finite(ce_df_inf)]
        df_cohort_agg <- if (length(valid_ce_df) > 0L) {
          as.integer(stats::median(valid_ce_df))
        } else {
          as.integer(df_fallback)
        }

        pvalue_cohort_agg <- 2 * stats::pt(
          abs(t_stat_cohort_agg), df = df_cohort_agg, lower.tail = FALSE
        )
        t_crit <- stats::qt(1 - alpha / 2, df = df_cohort_agg)
        ci_cohort_agg <- c(
          att_cohort_agg - t_crit * se_cohort_agg,
          att_cohort_agg + t_crit * se_cohort_agg
        )
      }
    }
  }

  # --- Step 5c: Top-level ATT fallback hierarchy ---
  att_top <- NA_real_
  se_top <- NA_real_

  if (aggregate == "overall" && !is.null(overall_effect)) {
    att_top <- overall_effect$att
    se_top <- overall_effect$se
  } else if (aggregate == "cohort" && !is.null(att_cohort_agg)) {
    att_top <- att_cohort_agg
    se_top <- se_cohort_agg
  } else {
    # aggregate="none" or "event_time": n_treated weighted average of (g,r) effects
    if (nrow(cohort_time_effects) > 0L &&
        any(!is.na(cohort_time_effects$att))) {
      valid_rows <- cohort_time_effects[
        !is.na(cohort_time_effects$att), , drop = FALSE
      ]
      w <- valid_rows$n_treated
      w_sum <- sum(w)
      if (w_sum > 0) {
        att_top <- stats::weighted.mean(valid_rows$att, w)
      }
    }
    se_top <- NA_real_  # Cannot derive overall SE from independent (g,r) SEs
  }

  # --- Step 5d: Degrees of freedom ---
  df_resid_summary <- NA_integer_
  df_inference_summary <- NA_integer_

  if (aggregate == "overall" && !is.null(overall_effect)) {
    df_resid_summary <- as.integer(overall_effect$df_resid)
    df_inference_summary <- as.integer(overall_effect$df_inference)
  } else if (aggregate == "cohort" && !is.null(cohort_effects)) {
    ce_df_inf <- vapply(cohort_effects, function(e) {
      as.numeric(e$df_inference)
    }, numeric(1))
    valid_df <- ce_df_inf[is.finite(ce_df_inf)]
    df_inference_summary <- if (length(valid_df) > 0L) {
      as.integer(stats::median(valid_df))
    } else {
      as.integer(df_fallback)
    }
    ce_df_res <- vapply(cohort_effects, function(e) {
      as.numeric(e$df_resid)
    }, numeric(1))
    valid_df_r <- ce_df_res[is.finite(ce_df_res)]
    df_resid_summary <- if (length(valid_df_r) > 0L) {
      as.integer(stats::median(valid_df_r))
    } else {
      as.integer(df_fallback)
    }
  } else {
    # none / event_time: median of (g,r) df
    valid_df_resid <- cohort_time_effects$df[
      is.finite(cohort_time_effects$df)
    ]
    df_resid_summary <- if (length(valid_df_resid) > 0L) {
      as.integer(stats::median(valid_df_resid))
    } else {
      as.integer(df_fallback)
    }
    valid_df_inf <- cohort_time_effects$df_inference[
      is.finite(cohort_time_effects$df_inference)
    ]
    df_inference_summary <- if (length(valid_df_inf) > 0L) {
      as.integer(stats::median(valid_df_inf))
    } else {
      as.integer(df_fallback)
    }
  }

  # --- n_control summary ---
  n_control_summary <- if (control_group_used == "never_treated") {
    as.integer(n_nt)
  } else {
    if (nrow(cohort_time_effects) > 0L) {
      as.integer(max(cohort_time_effects$n_control, na.rm = TRUE))
    } else {
      as.integer(n_nt)
    }
  }

  # --- n_treated summary ---
  n_treated_summary <- as.integer(validated$n_treated)

  # --- top-level K summary ---
  top_level_K <- NULL
  if (nrow(cohort_time_effects) > 0L &&
      "K" %in% names(cohort_time_effects)) {
    valid_k <- cohort_time_effects$K[is.finite(cohort_time_effects$K)]
    if (length(valid_k) > 0L) {
      top_level_K <- as.integer(max(valid_k))
    }
  }

  # --- Build lwdid_result ---
  result <- new_lwdid_result(
    att = att_top,
    se_att = se_top,
    t_stat = NA_real_,
    pvalue = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    df_resid = df_resid_summary,
    df_inference = df_inference_summary,
    nobs = validated$n_obs,
    n_treated = n_treated_summary,
    n_control = n_control_summary,
    K = top_level_K,
    tpost1 = NA_integer_,
    depvar = y,
    rolling = rolling,
    estimator = validated$validated_params$estimator %||% "ra",
    method = "staggered",
    vce_type = if (is.null(vce)) "homoskedastic" else vce,
    cluster_var = cluster_var,
    n_clusters = NULL,
    alpha = alpha,
    is_staggered = TRUE,
    controls_used = !is.null(controls),
    controls = controls %||% character(0),
    include_pretreatment =
      validated$validated_params$include_pretreatment,
    control_group = validated$validated_params$control_group,
    control_group_used = control_group_used,
    params = NULL, bse = NULL, vcov_matrix = NULL,
    resid = NULL, data = NULL,
    att_by_period = NULL, att_pre_treatment = NULL,
    parallel_trends_test = NULL,
    aggregate = aggregate,
    cohorts = valid_cohorts,
    cohort_sizes = cohort_sizes,
    n_never_treated = as.integer(n_nt),
    att_by_cohort = if (!is.null(cohort_effects) && length(cohort_effects) > 0L) {
      do.call(rbind, lapply(cohort_effects, function(ce) {
        data.frame(cohort = ce$cohort, att = ce$att, se = ce$se,
                   ci_lower = ce$ci_lower, ci_upper = ce$ci_upper,
                   pvalue = ce$pvalue, n_units = ce$n_units,
                   n_periods = ce$n_periods,
                   stringsAsFactors = FALSE)
      }))
    } else { NULL },
    att_by_cohort_time = cohort_time_effects,
    cohort_effects = cohort_effects,
    att_overall = if (!is.null(overall_effect)) overall_effect$att else NULL,
    se_overall = if (!is.null(overall_effect)) overall_effect$se else NULL,
    ci_overall_lower = if (!is.null(overall_effect)) overall_effect$ci_lower else NULL,
    ci_overall_upper = if (!is.null(overall_effect)) overall_effect$ci_upper else NULL,
    t_stat_overall = if (!is.null(overall_effect)) overall_effect$t_stat else NULL,
    pvalue_overall = if (!is.null(overall_effect)) overall_effect$pvalue else NULL,
    event_time_effects = event_time_effects,
    cohort_weights = NULL,
    ri_pvalue = NULL, ri_distribution = NULL,
    ri_seed = NULL, rireps = NULL, ri_method = NULL,
    ri_valid = NULL, ri_failed = NULL,
    ri_error = NULL, ri_target = NULL,
    diagnostics = NULL, warning_diagnostics = list(),
    propensity_scores = NULL,
    matched_data = NULL, n_matched = NULL,
    match_rate = NULL, weights_cv = NULL,
    warnings_log = list(),
    call = NULL, lwdid_version = NULL,
    ivar = ivar, tvar = tvar,
    is_quarterly = validated$is_quarterly %||% FALSE
  )

  # Store extra fields
  result$exclude_pre_periods <- exclude_pre_periods
  result$cohort_time_effects <- result$att_by_cohort_time
  result$overall_effect <- overall_effect
  result$att_cohort_agg <- att_cohort_agg
  result$se_cohort_agg <- se_cohort_agg
  result$t_stat_cohort_agg <- t_stat_cohort_agg
  result$pvalue_cohort_agg <- pvalue_cohort_agg
  if (!is.null(ci_cohort_agg)) {
    result$ci_lower_cohort_agg <- ci_cohort_agg[1]
    result$ci_upper_cohort_agg <- ci_cohort_agg[2]
  } else {
    result$ci_lower_cohort_agg <- NA_real_
    result$ci_upper_cohort_agg <- NA_real_
  }

  # Compute top-level inference if att and se are both finite
  top_level_inference_dist <- .resolve_top_level_inference_dist(
    validated$validated_params$estimator %||% "ra"
  )
  result$inference_dist <- top_level_inference_dist

  if (is.finite(att_top) && is.finite(se_top) && se_top > 0) {
    result$t_stat <- att_top / se_top

    if (identical(top_level_inference_dist, "normal")) {
      z_crit <- stats::qnorm(1 - alpha / 2)
      result$pvalue <- 2 * (1 - stats::pnorm(abs(result$t_stat)))
      result$ci_lower <- att_top - z_crit * se_top
      result$ci_upper <- att_top + z_crit * se_top
    } else if (is.finite(df_inference_summary) && df_inference_summary > 0L) {
      result$pvalue <- 2 * stats::pt(
        abs(result$t_stat), df = df_inference_summary, lower.tail = FALSE
      )
      t_crit <- stats::qt(1 - alpha / 2, df = df_inference_summary)
      result$ci_lower <- att_top - t_crit * se_top
      result$ci_upper <- att_top + t_crit * se_top
    }
  }

  result
}
