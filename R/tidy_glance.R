#' @title Tidy and Glance methods for lwdid_result (modelsummary / broom compatibility)
#' @description
#' S3 methods making lwdid_result objects compatible with the broom / modelsummary
#' ecosystem. Also provides fixest-style helper accessors.
#' @name tidy_glance
#' @family lwdid-results
NULL

#' Re-export tidy generic from generics
#' @param x An object to tidy
#' @param ... Additional arguments passed to methods
#' @importFrom generics tidy
#' @export
tidy <- generics::tidy

#' Re-export glance generic from generics
#' @param x An object to glance at
#' @param ... Additional arguments passed to methods
#' @importFrom generics glance
#' @export
glance <- generics::glance

# ===========================================================================
# tidy.lwdid_result
# ===========================================================================

#' Tidy a lwdid_result object
#'
#' @description
#' Extract estimation results in a tidy data.frame format compatible
#' with modelsummary and broom ecosystem. Follows Lee & Wooldridge (2025)
#' ATT definitions.
#'
#' @param x A lwdid_result object
#' @param type Character. Level of detail:
#'   \code{"overall"} (default) --- single row with overall ATT;
#'   \code{"effects"} --- all (g,r) pair effects from staggered estimation;
#'   \code{"event_time"} --- event-time aggregated effects.
#' @param conf.int Logical. Include confidence intervals? Default TRUE.
#' @param conf.level Numeric. Confidence level. Default 0.95.
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame with columns:
#'   \code{term}, \code{estimate}, \code{std.error}, \code{statistic},
#'   \code{p.value}, and optionally \code{conf.low}, \code{conf.high}.
#'
#' @seealso \code{\link{glance.lwdid_result}}, \code{\link{se}},
#'   \code{\link{compare}}
#' @family lwdid-results
#'
#' @examples
#' \donttest{
#'   # Create a minimal lwdid_result for demonstration
#'   res <- structure(list(
#'     att = 0.5, se_att = 0.1, t_stat = 5.0, pvalue = 0.001,
#'     ci_lower = 0.3, ci_upper = 0.7, nobs = 100L,
#'     n_treated = 50L, n_control = 50L, rolling = "demean",
#'     estimator = "ra", vce_type = "HC3", method = "proc21"
#'   ), class = "lwdid_result")
#'
#'   # Overall ATT in tidy format
#'   tidy(res)
#'
#'   # Model-level statistics
#'   glance(res)
#' }
#'
#' @importFrom generics tidy
#' @export
tidy.lwdid_result <- function(x, type = c("overall", "effects", "event_time"),
                              conf.int = TRUE, conf.level = 0.95, ...) {
  type <- match.arg(type)

  if (type == "overall") {
    out <- data.frame(
      term      = "ATT",
      estimate  = x$att,
      std.error = x$se_att,
      statistic = x$t_stat,
      p.value   = x$pvalue,
      stringsAsFactors = FALSE
    )
    if (conf.int) {
      out$conf.low  <- x$ci_lower
      out$conf.high <- x$ci_upper
    }

  } else if (type == "effects") {
    # Prefer att_by_cohort_time (data.frame with cohort, period, att, se, ...)
    eff <- x$att_by_cohort_time
    if (is.null(eff) || !is.data.frame(eff) || nrow(eff) == 0L) {
      return(data.frame(
        term = character(0), estimate = numeric(0),
        std.error = numeric(0), statistic = numeric(0),
        p.value = numeric(0), stringsAsFactors = FALSE
      ))
    }
    t_vals <- if (!is.null(eff$t_stat)) eff$t_stat else eff$att / eff$se
    p_vals <- if (!is.null(eff$pvalue)) eff$pvalue else {
      2 * stats::pnorm(-abs(t_vals))
    }
    out <- data.frame(
      term      = paste0("ATT(g=", eff$cohort, ",t=", eff$period, ")"),
      estimate  = eff$att,
      std.error = eff$se,
      statistic = t_vals,
      p.value   = p_vals,
      stringsAsFactors = FALSE
    )
    if (conf.int) {
      crit <- stats::qnorm(1 - (1 - conf.level) / 2)
      out$conf.low  <- if (!is.null(eff$ci_lower)) eff$ci_lower else eff$att - crit * eff$se
      out$conf.high <- if (!is.null(eff$ci_upper)) eff$ci_upper else eff$att + crit * eff$se
    }

  } else {
    # event_time
    et <- tryCatch(
      extract_effects(x, type = "event_time"),
      error = function(e) NULL
    )
    if (is.null(et) || !is.data.frame(et) || nrow(et) == 0L) {
      return(data.frame(
        term = character(0), estimate = numeric(0),
        std.error = numeric(0), statistic = numeric(0),
        p.value = numeric(0), stringsAsFactors = FALSE
      ))
    }
    t_vals <- if (!is.null(et$t_stat)) et$t_stat else et$att / et$se
    p_vals <- if (!is.null(et$pvalue)) et$pvalue else {
      2 * stats::pnorm(-abs(t_vals))
    }
    out <- data.frame(
      term      = paste0("event_time=", et$event_time),
      estimate  = et$att,
      std.error = et$se,
      statistic = t_vals,
      p.value   = p_vals,
      stringsAsFactors = FALSE
    )
    if (conf.int) {
      crit <- stats::qnorm(1 - (1 - conf.level) / 2)
      out$conf.low  <- if (!is.null(et$ci_lower)) et$ci_lower else et$att - crit * et$se
      out$conf.high <- if (!is.null(et$ci_upper)) et$ci_upper else et$att + crit * et$se
    }
  }

  out
}

# ===========================================================================
# glance.lwdid_result
# ===========================================================================

#' Glance at a lwdid_result object
#'
#' @description
#' Extract model-level statistics in a single-row data.frame compatible
#' with the broom / modelsummary ecosystem.
#'
#' @param x A lwdid_result object
#' @param ... Additional arguments (ignored).
#'
#' @return A single-row data.frame with columns:
#'   \code{nobs}, \code{n.treated}, \code{n.control}, \code{transformation},
#'   \code{estimator}, \code{vce.type}, \code{method}, \code{df.residual}.
#'
#' @seealso \code{\link{tidy.lwdid_result}}, \code{\link{se}},
#'   \code{\link{compare}}
#' @family lwdid-results
#'
#' @importFrom generics glance
#' @export
glance.lwdid_result <- function(x, ...) {
  data.frame(
    nobs           = if (!is.null(x$nobs)) x$nobs else NA_integer_,
    n.treated      = if (!is.null(x$n_treated)) x$n_treated else NA_integer_,
    n.control      = if (!is.null(x$n_control)) x$n_control else NA_integer_,
    transformation = if (!is.null(x$rolling)) x$rolling else NA_character_,
    estimator      = if (!is.null(x$estimator)) x$estimator else "ra",
    vce.type       = if (!is.null(x$vce_type)) x$vce_type else NA_character_,
    method         = if (!is.null(x$method)) x$method else NA_character_,
    df.residual    = if (!is.null(x$df_resid)) x$df_resid else NA_integer_,
    stringsAsFactors = FALSE
  )
}

# ===========================================================================
# fixest-style helper methods
# ===========================================================================

#' Extract standard errors
#'
#' @description
#' Generic function following fixest package conventions for extracting
#' standard errors from model objects. For lwdid_result objects, returns
#' the standard errors of the ATT estimate(s).
#'
#' @param x A model object.
#' @param ... Additional arguments passed to methods.
#' @return Named numeric vector of standard errors.
#'
#' @seealso \code{\link{tstat}}, \code{\link{pvalue}}, \code{\link{coeftable}},
#'   \code{\link{tidy.lwdid_result}}
#' @family lwdid-fixest
#' @export
se <- function(x, ...) UseMethod("se")

#' Extract t-statistics
#'
#' @description
#' Generic function following fixest package conventions for extracting
#' t-statistics from model objects. For lwdid_result objects, returns
#' the t-statistics of the ATT estimate(s).
#'
#' @param x A model object.
#' @param ... Additional arguments passed to methods.
#' @return Named numeric vector of t-statistics.
#'
#' @seealso \code{\link{se}}, \code{\link{pvalue}}, \code{\link{coeftable}},
#'   \code{\link{tidy.lwdid_result}}
#' @family lwdid-fixest
#' @export
tstat <- function(x, ...) UseMethod("tstat")

#' Extract p-values
#'
#' @description
#' Generic function following fixest package conventions for extracting
#' p-values from model objects. For lwdid_result objects, returns
#' the p-values of the ATT estimate(s).
#'
#' @param x A model object.
#' @param ... Additional arguments passed to methods.
#' @return Named numeric vector of p-values.
#'
#' @seealso \code{\link{se}}, \code{\link{tstat}}, \code{\link{coeftable}},
#'   \code{\link{tidy.lwdid_result}}
#' @family lwdid-fixest
#' @export
pvalue <- function(x, ...) UseMethod("pvalue")

#' Extract coefficient table
#'
#' @description
#' Generic function following fixest package conventions for extracting
#' a full coefficient table from model objects. For lwdid_result objects,
#' returns a matrix with estimates, standard errors, t-values, and p-values.
#'
#' @param x A model object.
#' @param ... Additional arguments passed to methods.
#' @return Matrix or data.frame with coefficient information.
#'
#' @seealso \code{\link{se}}, \code{\link{tstat}}, \code{\link{pvalue}},
#'   \code{\link{tidy.lwdid_result}}
#' @family lwdid-fixest
#' @export
coeftable <- function(x, ...) UseMethod("coeftable")

#' @rdname se
#' @examples
#' \donttest{
#'   res <- structure(list(
#'     att = 0.5, se_att = 0.1, t_stat = 5.0, pvalue = 0.001,
#'     ci_lower = 0.3, ci_upper = 0.7, nobs = 100L,
#'     n_treated = 50L, n_control = 50L, rolling = "demean",
#'     estimator = "ra", vce_type = "HC3", method = "proc21"
#'   ), class = "lwdid_result")
#'
#'   se(res)        # Standard errors
#'   tstat(res)     # t-statistics
#'   pvalue(res)    # p-values
#'   coeftable(res) # Full coefficient table
#' }
#' @export
se.lwdid_result <- function(x, ...) {
  if (!is.null(x$bse) && length(x$bse) > 0L) {
    return(x$bse)
  }
  stats::setNames(x$se_att, "ATT")
}

#' @rdname tstat
#' @export
tstat.lwdid_result <- function(x, ...) {
  if (!is.null(x$params) && length(x$params) > 0L &&
      !is.null(x$bse) && length(x$bse) > 0L) {
    return(x$params / x$bse)
  }
  stats::setNames(x$t_stat, "ATT")
}

#' @rdname pvalue
#' @export
pvalue.lwdid_result <- function(x, ...) {
  if (!is.null(x$params) && length(x$params) > 0L &&
      !is.null(x$bse) && length(x$bse) > 0L) {
    tvals <- x$params / x$bse
    df_val <- if (!is.null(x$df_resid) && is.finite(x$df_resid)) x$df_resid else Inf
    pvals <- 2 * stats::pt(-abs(tvals), df = df_val)
    names(pvals) <- names(x$params)
    return(pvals)
  }
  stats::setNames(x$pvalue, "ATT")
}

#' @rdname coeftable
#' @export
coeftable.lwdid_result <- function(x, ...) {
  est <- coef(x)
  se_val <- se.lwdid_result(x)

  # Ensure vectors are compatible

  if (length(se_val) != length(est)) {
    # Fallback to overall ATT
    est <- stats::setNames(x$att, "ATT")
    se_val <- stats::setNames(x$se_att, "ATT")
  }

  tval <- est / se_val
  df_val <- if (!is.null(x$df_resid) && is.finite(x$df_resid)) x$df_resid else Inf
  pval <- 2 * stats::pt(-abs(tval), df = df_val)

  out <- cbind(
    Estimate    = est,
    `Std. Error` = se_val,
    `t value`   = tval,
    `Pr(>|t|)`  = pval
  )
  rownames(out) <- names(est)
  out
}
