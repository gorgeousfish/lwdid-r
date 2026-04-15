pkg_root <- normalizePath("/Users/cxy/Desktop/lwdid_r/lwdid-r", winslash = "/")
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) == 0L) {
  stop("This audit must be run via Rscript --file.", call. = FALSE)
}
script_dir <- normalizePath(dirname(sub("^--file=", "", file_arg[[1L]])), winslash = "/")

failures <- character()

record_failure <- function(label, detail) {
  failures <<- c(failures, sprintf("[%s] %s", label, detail))
}

run_without_warnings <- function(expr, label) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      record_failure(label, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}

static_probe <- file.path(script_dir, "e8_04_package_hardening_static_regression.R")
tryCatch(
  source(static_probe, local = new.env(parent = globalenv())),
  error = function(e) record_failure("static-hardening", conditionMessage(e))
)

test_file <- file.path(pkg_root, "tests", "testthat", "test-sensitivity-pre-period.R")
tryCatch(
  invisible(parse(file = test_file)),
  error = function(e) record_failure("parse-test", conditionMessage(e))
)

rd_files <- c(
  file.path(pkg_root, "man", "aggregate_to_overall.Rd"),
  file.path(pkg_root, "man", "compute_inference.Rd")
)
for (rd_file in rd_files) {
  run_without_warnings(
    tryCatch(
      tools::parse_Rd(rd_file),
      error = function(e) record_failure("parse-rd", paste(basename(rd_file), conditionMessage(e)))
    ),
    label = paste("parse-rd", basename(rd_file))
  )
}

temp_lib <- tempfile("lwdid-hardening-lib-")
dir.create(temp_lib, recursive = TRUE)
build_dir <- tempfile("lwdid-hardening-build-")
dir.create(build_dir, recursive = TRUE)

install_cmd <- file.path(R.home("bin"), "R")
build_args <- c(
  "CMD", "build",
  "--no-manual",
  pkg_root
)
old_wd <- setwd(build_dir)

build_result <- system2(install_cmd, build_args, stdout = TRUE, stderr = TRUE)
build_status <- attr(build_result, "status")
if (!is.null(build_status) && build_status != 0L) {
  record_failure(
    "build",
    paste(build_result, collapse = "\n")
  )
} else {
  tarball <- Sys.glob(file.path(build_dir, "lwdid_*.tar.gz"))
  if (length(tarball) != 1L) {
    record_failure("build", "expected exactly one source tarball")
  } else {
    tar_contents <- utils::untar(tarball, list = TRUE)
    top_level_entries <- basename(tar_contents)

    leaked_tmp <- grep("^\\.tmp_", top_level_entries, value = TRUE)
    if (length(leaked_tmp) > 0L) {
      record_failure(
        "tmp-files",
        paste(sort(unique(leaked_tmp)), collapse = ", ")
      )
    }

    leaked_non_standard <- intersect(top_level_entries, c("_check_desc.R", "_manual_check.R"))
    if (length(leaked_non_standard) > 0L) {
      record_failure(
        "top-level-files",
        paste(sort(unique(leaked_non_standard)), collapse = ", ")
      )
    }
  }

  if (length(tarball) == 1L) {
    install_args <- c(
      "CMD", "INSTALL",
      "-l", temp_lib,
      tarball
    )
    install_result <- system2(install_cmd, install_args, stdout = TRUE, stderr = TRUE)
    install_status <- attr(install_result, "status")
    if (!is.null(install_status) && install_status != 0L) {
      record_failure(
        "install",
        paste(install_result, collapse = "\n")
      )
    } else {
      load_expr <- sprintf("loadNamespace(%s, lib.loc = %s)", dQuote("lwdid"), dQuote(temp_lib))
      load_cmd <- sprintf(
        "R_DEFAULT_PACKAGES=NULL %s --vanilla -q -e %s",
        shQuote(install_cmd),
        shQuote(load_expr)
      )
      load_result <- system2(
        "/bin/sh",
        c("-c", load_cmd),
        stdout = TRUE,
        stderr = TRUE
      )
      load_status <- attr(load_result, "status")
      if (!is.null(load_status) && load_status != 0L) {
        record_failure(
          "namespace-load",
          paste(load_result, collapse = "\n")
        )
      }
    }
  }
}

setwd(old_wd)
unlink(temp_lib, recursive = TRUE, force = TRUE)
unlink(build_dir, recursive = TRUE, force = TRUE)

if (length(failures) > 0L) {
  stop(
    paste(
      c("Package hardening regression failed:", failures),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

message("Package hardening regression passed.")
