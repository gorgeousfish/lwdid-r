args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1L) args[[1L]] else "/Users/cxy/Desktop/lwdid_r"
pkg_root <- file.path(repo_root, "lwdid-r")

namespace_path <- file.path(pkg_root, "NAMESPACE")
zzz_path <- file.path(pkg_root, "R", "zzz.R")

if (!file.exists(namespace_path)) {
  stop("NAMESPACE not found: ", namespace_path, call. = FALSE)
}
if (!file.exists(zzz_path)) {
  stop("zzz.R not found: ", zzz_path, call. = FALSE)
}

namespace_lines <- readLines(namespace_path, warn = FALSE, encoding = "UTF-8")
zzz_text <- paste(readLines(zzz_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

required_namespace_lines <- c(
  "importFrom(stats,as.formula)",
  "importFrom(stats,binomial)",
  "importFrom(stats,coef)",
  "importFrom(stats,complete.cases)",
  "importFrom(stats,glm)",
  "importFrom(stats,hatvalues)",
  "importFrom(stats,lm)",
  "importFrom(stats,lm.fit)",
  "importFrom(stats,model.frame)",
  "importFrom(stats,pnorm)",
  "importFrom(stats,predict)",
  "importFrom(stats,pt)",
  "importFrom(stats,qnorm)",
  "importFrom(stats,qt)",
  "importFrom(stats,residuals)",
  "importFrom(stats,rnorm)",
  "importFrom(stats,runif)",
  "importFrom(stats,sd)",
  "importFrom(stats,setNames)",
  "importFrom(stats,var)",
  "importFrom(stats,vcov)",
  "importFrom(utils,head)",
  "importFrom(utils,tail)"
)

required_globals <- c(
  "D_ig",
  "N",
  "alpha_hat",
  "beta_hat",
  "pre_mean",
  "ref_val",
  "y_trans",
  "y_trans_pre"
)

missing_namespace <- required_namespace_lines[!required_namespace_lines %in% namespace_lines]
missing_globals <- required_globals[
  !vapply(
    required_globals,
    function(symbol) grepl(sprintf('"%s"', symbol), zzz_text, fixed = TRUE),
    logical(1)
  )
]

issues <- character(0)
if (!grepl(".datatable.aware <- TRUE", zzz_text, fixed = TRUE)) {
  issues <- c(issues, "missing .datatable.aware <- TRUE sentinel in R/zzz.R")
}
if (length(missing_namespace) > 0L) {
  issues <- c(
    issues,
    paste0(
      "missing namespace imports: ",
      paste(missing_namespace, collapse = ", ")
    )
  )
}
if (length(missing_globals) > 0L) {
  issues <- c(
    issues,
    paste0(
      "missing globalVariables symbols: ",
      paste(missing_globals, collapse = ", ")
    )
  )
}

if (length(issues) > 0L) {
  writeLines(c("R code problems regression: FAIL", paste0("- ", issues)))
  quit(status = 1L)
}

writeLines("R code problems regression: PASS")
