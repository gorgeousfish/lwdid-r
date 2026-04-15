# control_groups.R
# Control group selection utilities for staggered DiD.
# Implements Never-Treated unit statistics for lw2025 Section 4.1
# and Procedure 4.1.
#
# Note: is_never_treated() and get_cohorts() are defined in utils.R.
# This file provides check_never_treated() for comprehensive NT
# statistics at the unit level.

# ============================================================================
# Never-Treated Statistics
# ============================================================================

#' @title Check Never-Treated unit presence and return statistics
#' @description Examines the panel data at the unit level to
#'   determine whether Never-Treated units exist and provides
#'   comprehensive statistics including NT count, treated count,
#'   and per-cohort sizes. Used by downstream functions to
#'   validate control group availability for aggregated effects
#'   (lw2026 equations 7.9-7.10 for cohort, 7.18-7.19 for overall).
#'
#' @param data data.table containing panel data
#' @param ivar character, name of unit identifier column
#' @param gvar character, name of cohort variable column
#' @return list with elements:
#'   \describe{
#'     \item{has_nt}{logical, whether any NT units exist}
#'     \item{n_nt}{integer, number of NT units}
#'     \item{n_treated}{integer, number of treated units}
#'     \item{cohort_sizes}{table, unit counts per cohort}
#'   }
#' @keywords internal
check_never_treated <- function(data, ivar, gvar) {
  # Extract unit-level cohort values (deduplicate to unit level)
  unit_g <- unique(data[, c(ivar, gvar), with = FALSE])
  g_vals <- unit_g[[gvar]]

  nt_mask <- is_never_treated(g_vals)
  cohort_vals <- g_vals[!nt_mask]

  list(
    has_nt = any(nt_mask),
    n_nt = sum(nt_mask),
    n_treated = sum(!nt_mask),
    cohort_sizes = table(cohort_vals)
  )
}

# ============================================================================
# Control Group Strategy Resolution
# ============================================================================

#' @title Resolve control group strategy
#' @description Resolves the control group strategy based on the
#'   combination of user-specified control_group, aggregate level,
#'   and presence of Never-Treated units. Implements the 8-case
#'   auto resolution logic and automatic switching for aggregation
#'   compatibility.
#'
#'   Auto resolution rules (lw2025 Section 4.1, Procedure 4.1):
#'   \itemize{
#'     \item has_nt + cohort/overall: "never_treated" (required for
#'       aggregated effects, lw2025 equations 7.9-7.10, 7.18-7.19)
#'     \item has_nt + none/event_time: "not_yet_treated" (event_time
#'       re-weights existing (g,r) effects, no new regressions needed,
#'       so NT-only control is not required)
#'     \item !has_nt + cohort/overall: error (lwdid_no_never_treated)
#'     \item !has_nt + none/event_time: "not_yet_treated"
#'   }
#'
#'   Automatic switching: when aggregate is "cohort" or "overall"
#'   and control_group is "not_yet_treated" or "all_others",
#'   switches to "never_treated" with a message notification.
#'
#'   Differences from Python reference implementation:
#'   \itemize{
#'     \item Returns list(resolved, switched) instead of plain string,
#'       adding the switched field for caller diagnostics.
#'     \item R supports event_time aggregate level. Python defines
#'       \code{VALID_AGGREGATE = ('none', 'cohort', 'overall')} and
#'       does not support event_time; R's event_time is a feature
#'       extension.
#'     \item Python raises SmallSampleWarning during control group
#'       resolution when NT count is low. R defers this check to
#'       the inference stage (VCE/bootstrap), where sample size
#'       diagnostics are more naturally handled.
#'   }
#'
#' @param control_group character, one of \code{"auto"},
#'   \code{"not_yet_treated"}, \code{"never_treated"},
#'   \code{"all_others"}
#' @param aggregate character, one of \code{"none"}, \code{"cohort"},
#'   \code{"overall"}, \code{"event_time"}
#' @param has_nt logical, whether Never-Treated units exist in the data
#' @return list with elements:
#'   \describe{
#'     \item{resolved}{character, the resolved control group strategy}
#'     \item{switched}{logical, whether automatic switching occurred
#'       (TRUE when the original strategy was overridden for
#'       aggregation compatibility)}
#'   }
#' @keywords internal
resolve_control_group <- function(control_group, aggregate, has_nt) {
  switched <- FALSE

  # Case 0: auto strategy resolution (8 combinations)
  if (control_group == "auto") {
    if (has_nt) {
      # Rules 1-2: NT units available
      if (aggregate %in% c("cohort", "overall")) {
        # Aggregated effects (eqs 7.9-7.10, 7.18-7.19) require NT controls
        return(list(resolved = "never_treated", switched = FALSE))
      } else {
        # aggregate %in% c("none", "event_time")
        # event_time re-weights existing (g,r) effects; no new regressions,
        # so NT-only control is not required. Python VALID_AGGREGATE does
        # not include event_time; R's support is a feature extension.
        return(list(resolved = "not_yet_treated", switched = FALSE))
      }
    } else {
      # Rules 3-4: AET (no NT units)
      if (aggregate %in% c("cohort", "overall")) {
        stop_lwdid(
          sprintf(
            paste("Cannot use '%s' aggregation without Never-Treated units.",
                  "Cohort/overall aggregation requires control_group = 'never_treated',",
                  "but no NT units found in data."),
            aggregate
          ),
          class = "lwdid_no_never_treated",
          aggregate = aggregate,
          control_group = "auto"
        )
      } else {
        return(list(resolved = "not_yet_treated", switched = FALSE))
      }
    }
  }

  # Case 1b: explicit never_treated but no NT units → error
  if (control_group == "never_treated" && !has_nt) {
    stop_lwdid(
      paste("control_group = 'never_treated' specified but no",
            "Never-Treated units found in data."),
      class = "lwdid_no_never_treated",
      aggregate = aggregate,
      control_group = "never_treated"
    )
  }

  # Case 2: automatic switching for aggregation compatibility
  if (aggregate %in% c("cohort", "overall") &&
      control_group %in% c("not_yet_treated", "all_others")) {
    # Check NT availability before switching
    if (!has_nt) {
      stop_lwdid(
        sprintf(
          paste("Cannot switch to 'never_treated' for '%s' aggregation:",
                "no Never-Treated units found in data."),
          aggregate
        ),
        class = "lwdid_no_never_treated",
        aggregate = aggregate,
        control_group = control_group
      )
    }
    message(sprintf(
      "[lwdid] Switching control_group from '%s' to 'never_treated' for '%s' aggregation.",
      control_group, aggregate
    ))
    return(list(resolved = "never_treated", switched = TRUE))
  }

  # Case 3: all_others warning for non-aggregation scenarios
  if (control_group == "all_others" && aggregate %in% c("none", "event_time")) {
    warn_lwdid(
      paste("control_group = 'all_others' includes already-treated units",
            "as controls, which may introduce bias. Consider using",
            "'not_yet_treated' instead."),
      class = "lwdid_data",
      detail = "all_others_bias"
    )
  }

  # Default: return as-is
  list(resolved = control_group, switched = switched)
}

# ============================================================================
# Control Group Mask Generation
# ============================================================================

#' @title Generate control group mask for a (g, r) cell
#' @description Generates a logical mask identifying valid control
#'   units for a specific (g, r) cell. Operates at the observation
#'   level on \code{data[[gvar]]} (does NOT use \code{unique()}).
#'
#'   \strong{FATAL-001}: \code{not_yet_treated} uses STRICT inequality
#'   \eqn{G_i > r}. This is critical for correct identification of
#'   control units per Procedure 4.1 Step 2 and equation 4.12 in
#'   lw2025. Using \eqn{\ge}{>=} would incorrectly include units
#'   treated at period r.
#'
#'   Theoretical basis — Theorem 4.1 derivation chain (lw2025):
#'   \itemize{
#'     \item Equation 4.7: defines the 2x2 DiD estimand for cohort g
#'       at reference period r, requiring untreated controls at both
#'       periods r and r+1.
#'     \item Equation 4.8: parallel trends assumption conditions on
#'       \eqn{G_i > r}, ensuring control units are not yet treated.
#'     \item Equation 4.9: no-anticipation assumption further requires
#'       that units with \eqn{G_i > r} have not altered behavior
#'       before their treatment date.
#'   }
#'   Together, equations 4.7-4.9 establish that only units with
#'   \eqn{G_i > r} (strict) are valid not-yet-treated controls.
#'   Units with \eqn{G_i = r} are treated starting at period r and
#'   violate the parallel trends requirement.
#'
#'   Calling convention: caller passes a single-period cross-section
#'   (one row per unit). Confirmed by Story E4-04: the (g,r) effect
#'   estimation loop first extracts \code{period_data}, then calls
#'   this function on \code{period_data}.
#'
#'   This is an internal function and does NOT perform basic input
#'   validation (empty data, column existence, gvar type). The
#'   caller is responsible for validation before calling.
#'
#' @param data data.table, single-period cross-section (one row per unit)
#' @param gvar character, name of cohort variable column in \code{data}
#' @param g integer, focal cohort value (units with \code{gvar == g}
#'   are the treated group, not controls)
#' @param r integer, reference period for the 2x2 DiD comparison
#' @param control_group character, one of \code{"not_yet_treated"},
#'   \code{"never_treated"}, \code{"all_others"}. Must NOT be
#'   \code{"auto"} — resolve via \code{resolve_control_group()} first.
#' @return logical vector of length \code{nrow(data)}, where
#'   \code{TRUE} indicates a valid control unit for this (g, r) cell
#' @keywords internal
get_valid_controls <- function(data, gvar, g, r, control_group) {
  g_vals <- data[[gvar]]

  # NT mask via is_never_treated() (handles NA, NaN, 0, Inf, near-zero)
  nt_mask <- is_never_treated(g_vals)

  if (control_group == "not_yet_treated") {
    # FATAL-001: STRICT inequality G_i > r
    # Procedure 4.1 Step 2, equation 4.12 (lw2025)
    # Units with G_i = r are treated at period r and MUST be excluded.
    # Theorem 4.1 derivation chain (eqs 4.7-4.9):
    #   eq 4.7: 2x2 DiD estimand requires untreated controls at r and r+1
    #   eq 4.8: parallel trends conditions on G_i > r (strict)
    #   eq 4.9: no-anticipation requires G_i > r units behave as untreated
    nyt_mask <- g_vals > r
    # NA safety: comparison with NA yields NA, coerce to FALSE
    nyt_mask[is.na(nyt_mask)] <- FALSE
    # Note: gvar=0 NT units: 0 > r is FALSE, but nt_mask captures them.
    # gvar=Inf NT units: Inf > r is TRUE, but nt_mask also captures them.
    # The union ensures correct handling regardless.
    mask <- nt_mask | nyt_mask

  } else if (control_group == "never_treated") {
    # Only NT units serve as controls
    mask <- nt_mask

  } else if (control_group == "all_others") {
    # All units except the focal cohort g
    # Explicit NA handling: NT units have NA/0/Inf gvar values,
    # g_vals != g may yield NA for NA gvar. Handle explicitly:
    # - NT units: always included (they are not cohort g)
    # - Non-NT units: included if g_vals != g
    #
    # Equivalence analysis: this could be simplified to
    #   is_never_treated(g_vals) | g_vals != g
    # but the explicit form is retained for three reasons:
    # 1. Readability: makes the NT special-case handling visible
    # 2. Defensive: avoids relying on NA comparison behavior
    # 3. Python alignment: matches the reference implementation structure
    mask <- (!nt_mask & g_vals != g) | nt_mask

  } else {
    # Invalid strategy (including "auto" which should be resolved first)
    stop_lwdid(
      sprintf(
        paste("Invalid control_group '%s' in get_valid_controls().",
              "Must be 'not_yet_treated', 'never_treated', or 'all_others'.",
              "Use resolve_control_group() to resolve 'auto' first."),
        control_group
      ),
      class = "lwdid_invalid_parameter",
      param = "control_group",
      value = control_group,
      allowed = c("not_yet_treated", "never_treated", "all_others")
    )
  }

  # Final NA safety net
  mask[is.na(mask)] <- FALSE

  # Empty control group check
  if (sum(mask) == 0L) {
    stop_lwdid(
      sprintf(
        "No valid control units for cohort g=%d, period r=%d with strategy '%s'.",
        g, r, control_group
      ),
      class = "lwdid_no_control",
      g = g, r = r, control_group = control_group
    )
  }

  mask
}
