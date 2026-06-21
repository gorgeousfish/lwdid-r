pkg_root <- normalizePath("/Users/cxy/Desktop/lwdid_r/lwdid-r", winslash = "/")
rd_file <- file.path(pkg_root, "man", "compute_inference.Rd")

failures <- character()

lines <- readLines(rd_file, warn = FALSE, encoding = "UTF-8")
if (any(grepl("\\\\%", lines, fixed = TRUE))) {
  failures <- c(
    failures,
    "Found double-escaped percent sequence '\\\\%' in compute_inference.Rd."
  )
}

conversion_warnings <- character()
conversion_result <- withCallingHandlers(
  tryCatch(
    capture.output(tools::Rd2txt(rd_file)),
    error = function(e) e
  ),
  warning = function(w) {
    conversion_warnings <<- c(conversion_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

if (inherits(conversion_result, "error")) {
  failures <- c(failures, conditionMessage(conversion_result))
}
if (length(conversion_warnings) > 0L) {
  failures <- c(failures, conversion_warnings)
}

if (length(failures) > 0L) {
  stop(
    paste(
      c("compute_inference.Rd regression failed:", failures),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

message("compute_inference.Rd regression passed.")
