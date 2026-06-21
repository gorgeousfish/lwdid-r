#' @title 比较多个 lwdid 估计结果
#' @description
#' 比较不同 specification（估计器、转换方法、控制组策略）下的 lwdid 估计结果，
#' 以标准化表格形式输出。基于 Lee & Wooldridge (2025, 2026) 的 ATT 定义。
#'
#' 不同 rolling 方法给出的是同一目标参数 tau 在不同识别假设下的估计，
#' compare() 的数学意义在于展示 ATT 在不同假设和方法下的稳健性。
#'
#' @param ... 命名的 lwdid_result 对象，或一个包含它们的列表
#' @param type 字符型。比较类型："overall"（默认，整体 ATT）或
#'   "effects"（分期效应）
#' @param stats 字符向量。包含的统计量。默认：
#'   c("estimate", "std.error", "ci", "p.value")
#' @param stars 逻辑值。是否添加显著性星号？默认 TRUE
#' @param digits 整数。小数位数。默认 3
#'
#' @return 一个 \code{lwdid_comparison} S3 对象（继承 data.frame），包含：
#'   \itemize{
#'     \item 每个模型 specification 为一列
#'     \item 统计量（估计值、标准误、CI、p 值等）为行
#'     \item 格式化的 print 方法
#'   }
#'
#' @seealso \code{\link{tidy.lwdid_result}}, \code{\link{print.lwdid_comparison}}
#' @family lwdid-results
#'
#' @export
#' @examples
#' \donttest{
#'   res1 <- structure(list(
#'     att = 0.5, se_att = 0.1, t_stat = 5.0, pvalue = 0.001,
#'     ci_lower = 0.3, ci_upper = 0.7, nobs = 100L,
#'     n_treated = 50L, n_control = 50L, rolling = "demean",
#'     estimator = "ra", vce_type = "HC3", method = "proc21"
#'   ), class = "lwdid_result")
#'
#'   res2 <- structure(list(
#'     att = 0.6, se_att = 0.15, t_stat = 4.0, pvalue = 0.005,
#'     ci_lower = 0.31, ci_upper = 0.89, nobs = 100L,
#'     n_treated = 50L, n_control = 50L, rolling = "demean",
#'     estimator = "ipw", vce_type = "HC3", method = "proc21"
#'   ), class = "lwdid_result")
#'
#'   # Compare two specifications
#'   compare(RA = res1, IPW = res2)
#'
#'   # Compare with effects breakdown
#'   compare(RA = res1, IPW = res2, type = "effects")
#' }
compare <- function(..., type = c("overall", "effects"),
                    stats = c("estimate", "std.error", "ci", "p.value"),
                    stars = TRUE, digits = 3L) {

  # --- 0. 参数匹配 ---
  type <- match.arg(type)
  digits <- as.integer(digits)
  stopifnot(is.logical(stars), length(stars) == 1L)
  stopifnot(is.numeric(digits), length(digits) == 1L, digits >= 0L)

  # --- 1. 收集模型对象 ---
  dots <- list(...)

  # 支持传入单个列表的情况

  if (length(dots) == 1L && is.list(dots[[1L]]) &&
      !inherits(dots[[1L]], "lwdid_result")) {
    dots <- dots[[1L]]
  }

  n_models <- length(dots)
  if (n_models == 0L) {
    stop("\u81f3\u5c11\u9700\u8981\u4e00\u4e2a lwdid_result \u5bf9\u8c61\u8fdb\u884c\u6bd4\u8f83", call. = FALSE)
  }

  # --- 2. 输入验证：确保所有对象都是 lwdid_result ---
  for (i in seq_len(n_models)) {
    if (!inherits(dots[[i]], "lwdid_result")) {
      stop(sprintf("\u7b2c %d \u4e2a\u53c2\u6570\u4e0d\u662f lwdid_result \u5bf9\u8c61", i),
           call. = FALSE)
    }
  }

  # --- 3. 自动命名 ---
  model_names <- names(dots)
  if (is.null(model_names) || any(!nzchar(model_names))) {
    # 为未命名模型自动生成描述性名称
    auto_names <- vapply(dots, function(m) {
      rolling_part <- if (!is.null(m$rolling) && !is.na(m$rolling)) {
        as.character(m$rolling)
      } else {
        "unknown"
      }
      est_part <- if (!is.null(m$estimator) && !is.na(m$estimator)) {
        as.character(m$estimator)
      } else {
        "est"
      }
      paste0(rolling_part, "_", est_part)
    }, character(1))

    if (is.null(model_names)) {
      model_names <- auto_names
    } else {
      # 仅填充空名称
      empty_idx <- !nzchar(model_names)
      model_names[empty_idx] <- auto_names[empty_idx]
    }

    # 处理重复名称
    if (anyDuplicated(model_names)) {
      model_names <- make.unique(model_names, sep = "_")
    }
  }

  # --- 4. 根据 type 提取统计量 ---
  if (type == "overall") {
    result <- .compare_overall(dots, model_names, stats, stars, digits)
  } else {
    result <- .compare_effects(dots, model_names, stats, stars, digits)
  }

  result
}


# =============================================================================
# 内部函数：比较 overall ATT
# =============================================================================
.compare_overall <- function(models, model_names, stats, stars, digits) {
  n_models <- length(models)

  # 提取各模型的核心统计量
  att_vals <- vapply(models, function(m) {
    if (!is.null(m$att) && is.finite(m$att)) m$att else NA_real_
  }, numeric(1))

  se_vals <- vapply(models, function(m) {
    if (!is.null(m$se_att) && is.finite(m$se_att)) m$se_att else NA_real_
  }, numeric(1))

  pval_vals <- vapply(models, function(m) {
    if (!is.null(m$pvalue) && is.finite(m$pvalue)) m$pvalue else NA_real_
  }, numeric(1))

  ci_lower_vals <- vapply(models, function(m) {
    if (!is.null(m$ci_lower) && is.finite(m$ci_lower)) m$ci_lower else NA_real_
  }, numeric(1))

  ci_upper_vals <- vapply(models, function(m) {
    if (!is.null(m$ci_upper) && is.finite(m$ci_upper)) m$ci_upper else NA_real_
  }, numeric(1))

  # 生成显著性星号
  star_strs <- if (stars) {
    vapply(pval_vals, function(p) {
      if (is.na(p)) return("")
      if (p < 0.01) return("***")
      if (p < 0.05) return("**")
      if (p < 0.10) return("*")
      ""
    }, character(1))
  } else {
    rep("", n_models)
  }

  # 构建输出表格行
  rows <- list()
  row_labels <- character(0)

  if ("estimate" %in% stats) {
    # ATT 估计值行（含星号）
    est_row <- vapply(seq_len(n_models), function(i) {
      if (is.na(att_vals[i])) return("NA")
      paste0(formatC(att_vals[i], format = "f", digits = digits), star_strs[i])
    }, character(1))
    rows <- c(rows, list(est_row))
    row_labels <- c(row_labels, "ATT")
  }

  if ("std.error" %in% stats) {
    # 标准误行（括号格式）
    se_row <- vapply(se_vals, function(s) {
      if (is.na(s)) return("")
      sprintf("(%s)", formatC(s, format = "f", digits = digits))
    }, character(1))
    rows <- c(rows, list(se_row))
    row_labels <- c(row_labels, "")
  }

  if ("ci" %in% stats) {
    # 置信区间行
    alpha_val <- models[[1]]$alpha %||% 0.05
    ci_level_pct <- round((1 - alpha_val) * 100)
    ci_row <- vapply(seq_len(n_models), function(i) {
      if (is.na(ci_lower_vals[i]) || is.na(ci_upper_vals[i])) return("")
      sprintf("[%s,%s]",
              formatC(ci_lower_vals[i], format = "f", digits = digits),
              formatC(ci_upper_vals[i], format = "f", digits = digits))
    }, character(1))
    rows <- c(rows, list(ci_row))
    row_labels <- c(row_labels, sprintf("%d%% CI", ci_level_pct))
  }

  if ("p.value" %in% stats) {
    # p 值行
    pval_row <- vapply(pval_vals, function(p) {
      if (is.na(p)) return("")
      if (p < 0.001) return("<0.001")
      formatC(p, format = "f", digits = digits)
    }, character(1))
    rows <- c(rows, list(pval_row))
    row_labels <- c(row_labels, "p-value")
  }

  # 模型信息行
  info_rows <- list()
  info_labels <- character(0)

  # N（观测数）
  n_row <- vapply(models, function(m) {
    n <- m$nobs
    if (!is.null(n) && is.finite(n)) as.character(as.integer(n)) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(n_row))
  info_labels <- c(info_labels, "N")

  # N_treated
  nt_row <- vapply(models, function(m) {
    n <- m$n_treated
    if (!is.null(n) && is.finite(n)) as.character(as.integer(n)) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(nt_row))
  info_labels <- c(info_labels, "N_treated")

  # N_control
  nc_row <- vapply(models, function(m) {
    n <- m$n_control
    if (!is.null(n) && is.finite(n)) as.character(as.integer(n)) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(nc_row))
  info_labels <- c(info_labels, "N_control")

  # Transformation（转换方法）
  trans_row <- vapply(models, function(m) {
    if (!is.null(m$rolling) && !is.na(m$rolling)) as.character(m$rolling) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(trans_row))
  info_labels <- c(info_labels, "Transformation")

  # Estimator（估计器）
  est_info_row <- vapply(models, function(m) {
    if (!is.null(m$estimator) && !is.na(m$estimator)) as.character(m$estimator) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(est_info_row))
  info_labels <- c(info_labels, "Estimator")

  # VCE
  vce_row <- vapply(models, function(m) {
    if (!is.null(m$vce_type) && !is.na(m$vce_type)) as.character(m$vce_type) else "NA"
  }, character(1))
  info_rows <- c(info_rows, list(vce_row))
  info_labels <- c(info_labels, "VCE")

  # Control group（控制组策略，仅在非全部相同时显示）
  cg_vals <- vapply(models, function(m) {
    if (!is.null(m$control_group) && !is.na(m$control_group)) {
      as.character(m$control_group)
    } else {
      "NA"
    }
  }, character(1))
  if (length(unique(cg_vals)) > 1L) {
    info_rows <- c(info_rows, list(cg_vals))
    info_labels <- c(info_labels, "Control group")
  }

  # 组装结果 data.frame
  all_rows <- c(rows, info_rows)
  all_labels <- c(row_labels, info_labels)

  # 构建 data.frame
  result_df <- as.data.frame(
    matrix(unlist(all_rows), nrow = length(all_rows), byrow = TRUE),
    stringsAsFactors = FALSE
  )
  names(result_df) <- model_names
  result_df <- cbind(data.frame(stat = all_labels, stringsAsFactors = FALSE),
                     result_df)
  rownames(result_df) <- NULL

  # 附加元数据
  attr(result_df, "n_models") <- n_models
  attr(result_df, "model_names") <- model_names
  attr(result_df, "type") <- "overall"
  attr(result_df, "stars") <- stars
  attr(result_df, "digits") <- digits
  attr(result_df, "n_stat_rows") <- length(rows)

  class(result_df) <- c("lwdid_comparison", "data.frame")
  result_df
}


# =============================================================================
# 内部函数：比较 effects（分期效应）
# =============================================================================
.compare_effects <- function(models, model_names, stats, stars, digits) {
  n_models <- length(models)

  # 尝试从 effects 字段提取分期效应

  effects_list <- lapply(models, function(m) {
    eff <- m$effects
    if (is.null(eff) || !is.data.frame(eff) || nrow(eff) == 0L) {
      # 尝试 cohort_effects 或 event_time_effects
      eff <- m$cohort_effects
      if (is.null(eff) || !is.data.frame(eff) || nrow(eff) == 0L) {
        eff <- m$event_time_effects
      }
    }
    eff
  })

  # 检查是否至少有一个模型有分期效应
  has_effects <- vapply(effects_list, function(e) {
    !is.null(e) && is.data.frame(e) && nrow(e) > 0L
  }, logical(1))

  if (!any(has_effects)) {
    stop("\u6ca1\u6709\u6a21\u578b\u5305\u542b\u5206\u671f\u6548\u5e94\uff08effects\uff09\u4fe1\u606f\uff0c\u65e0\u6cd5\u8fdb\u884c effects \u6bd4\u8f83",
         call. = FALSE)
  }

  # 取第一个有效应的模型来确定行（g,r 对或 event_time）
  ref_idx <- which(has_effects)[1L]
  ref_eff <- effects_list[[ref_idx]]

  # 尝试确定行标识列
  id_cols <- intersect(names(ref_eff), c("g", "r", "cohort", "time", "event_time"))
  if (length(id_cols) == 0L) {
    # 使用行号作为标识
    id_cols <- NULL
    period_labels <- paste0("Period_", seq_len(nrow(ref_eff)))
  } else {
    period_labels <- apply(ref_eff[, id_cols, drop = FALSE], 1, function(row) {
      paste(row, collapse = ",")
    })
  }

  n_periods <- length(period_labels)

  # 对每个模型/每个期间提取 att 和 se
  result_rows <- list()
  result_labels <- character(0)

  for (p in seq_len(n_periods)) {
    # ATT 行
    att_row <- vapply(seq_len(n_models), function(i) {
      eff <- effects_list[[i]]
      if (is.null(eff) || nrow(eff) < p) return("NA")
      att_val <- eff$att[p]
      pval <- if ("pvalue" %in% names(eff)) eff$pvalue[p] else NA_real_
      star_str <- ""
      if (stars && !is.na(pval)) {
        if (pval < 0.01) star_str <- "***"
        else if (pval < 0.05) star_str <- "**"
        else if (pval < 0.10) star_str <- "*"
      }
      if (is.na(att_val)) return("NA")
      paste0(formatC(att_val, format = "f", digits = digits), star_str)
    }, character(1))
    result_rows <- c(result_rows, list(att_row))
    result_labels <- c(result_labels, period_labels[p])

    # SE 行（如果需要）
    if ("std.error" %in% stats) {
      se_row <- vapply(seq_len(n_models), function(i) {
        eff <- effects_list[[i]]
        if (is.null(eff) || nrow(eff) < p) return("")
        se_val <- if ("se" %in% names(eff)) eff$se[p] else NA_real_
        if (is.na(se_val)) return("")
        sprintf("(%s)", formatC(se_val, format = "f", digits = digits))
      }, character(1))
      result_rows <- c(result_rows, list(se_row))
      result_labels <- c(result_labels, "")
    }
  }

  # 组装 data.frame
  result_df <- as.data.frame(
    matrix(unlist(result_rows), nrow = length(result_rows), byrow = TRUE),
    stringsAsFactors = FALSE
  )
  names(result_df) <- model_names
  result_df <- cbind(data.frame(stat = result_labels, stringsAsFactors = FALSE),
                     result_df)
  rownames(result_df) <- NULL

  # 附加元数据
  attr(result_df, "n_models") <- n_models
  attr(result_df, "model_names") <- model_names
  attr(result_df, "type") <- "effects"
  attr(result_df, "stars") <- stars
  attr(result_df, "digits") <- digits
  attr(result_df, "n_periods") <- n_periods

  class(result_df) <- c("lwdid_comparison", "data.frame")
  result_df
}


# =============================================================================
# S3 方法：格式化打印
# =============================================================================

#' Print method for lwdid_comparison objects
#'
#' @description
#' Formats and prints a comparison table showing side-by-side results
#' from multiple lwdid specifications.
#'
#' @param x An lwdid_comparison object returned by \code{\link{compare}}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#'
#' @seealso \code{\link{compare}}
#' @family lwdid-results
#' @export
print.lwdid_comparison <- function(x, ...) {
  type <- attr(x, "type") %||% "overall"
  model_names <- attr(x, "model_names") %||% names(x)[-1L]
  stars_used <- attr(x, "stars") %||% TRUE
  n_stat_rows <- attr(x, "n_stat_rows") %||% 0L

  # 计算各列宽度
  n_models <- length(model_names)
  col_widths <- vapply(seq_len(n_models), function(i) {
    col_data <- x[[i + 1L]]  # 跳过 stat 列
    max(nchar(model_names[i]), max(nchar(as.character(col_data)), na.rm = TRUE))
  }, integer(1))
  # 确保最小宽度

  col_widths <- pmax(col_widths, 10L)

  label_width <- max(nchar(x$stat), 16L)

  # 总表格宽度
  total_width <- label_width + 2L + sum(col_widths + 2L)

  # 标题
  cat("\nlwdid Comparison Table\n")
  cat(strrep("=", total_width), "\n")

  # 表头行
  header <- sprintf("%-*s", label_width, "")
  for (i in seq_len(n_models)) {
    header <- paste0(header, "  ", sprintf("%-*s", col_widths[i], model_names[i]))
  }
  cat(header, "\n")
  cat(strrep("-", total_width), "\n")

  # 统计量行
  n_rows <- nrow(x)
  for (r in seq_len(n_rows)) {
    row_label <- x$stat[r]

    # 在统计量行和信息行之间添加分隔线
    if (n_stat_rows > 0L && r == n_stat_rows + 1L) {
      cat(strrep("-", total_width), "\n")
    }

    line <- sprintf("%-*s", label_width, row_label)
    for (i in seq_len(n_models)) {
      cell <- as.character(x[[i + 1L]][r])
      line <- paste0(line, "  ", sprintf("%-*s", col_widths[i], cell))
    }
    cat(line, "\n")
  }

  # 底部
  cat(strrep("=", total_width), "\n")
  if (stars_used) {
    cat("*** p<0.01, ** p<0.05, * p<0.1\n")
  }
  cat("\n")

  invisible(x)
}
