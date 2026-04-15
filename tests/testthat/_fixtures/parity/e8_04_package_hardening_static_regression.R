pkg_root <- normalizePath("/Users/cxy/Desktop/lwdid_r/lwdid-r", winslash = "/")

failures <- character()

record_failure <- function(label, detail) {
  failures <<- c(failures, sprintf("[%s] %s", label, detail))
}

read_description_fields <- function(path, field) {
  dcf <- read.dcf(path)
  value <- dcf[1, field]
  if (is.na(value) || !nzchar(value)) {
    return(character())
  }

  entries <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  gsub("\\s*\\(.*\\)$", "", entries)
}

extract_namespace_packages <- function(paths) {
  packages <- character()

  for (path in paths) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lines <- sub("#.*$", "", lines)
    lines <- lines[nzchar(trimws(lines))]
    matches <- gregexpr("([A-Za-z][A-Za-z0-9.]*):::{0,1}", lines, perl = TRUE)
    found <- regmatches(lines, matches)
    if (length(found) == 0L) {
      next
    }

    tokens <- unlist(found, use.names = FALSE)
    if (length(tokens) == 0L) {
      next
    }

    tokens <- sub(":::{0,1}$", "", tokens)
    tokens <- tokens[!grepl("\\.R$", tokens, perl = TRUE)]
    packages <- c(packages, tokens)
  }

  unique(packages)
}

find_non_ascii_files <- function(paths) {
  offenders <- character()

  for (path in paths) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    code_lines <- sub("#.*$", "", lines)
    code_lines <- code_lines[nzchar(trimws(code_lines))]
    if (any(grepl("[^\u0001-\u007F]", code_lines, perl = TRUE))) {
      offenders <- c(offenders, basename(path))
    }
  }

  unique(offenders)
}

description_path <- file.path(pkg_root, "DESCRIPTION")
pkg_name <- read.dcf(description_path)[1, "Package"]
imports <- read_description_fields(description_path, "Imports")
depends <- read_description_fields(description_path, "Depends")
suggests <- read_description_fields(description_path, "Suggests")

r_declared <- c(imports, depends, suggests, "base")
r_files <- Sys.glob(file.path(pkg_root, "R", "*.R"))
test_files <- Sys.glob(file.path(pkg_root, "tests", "testthat", "*.R"))

non_ascii <- find_non_ascii_files(r_files)
if (length(non_ascii) > 0L) {
  record_failure(
    "non-ascii-r",
    paste(sort(non_ascii), collapse = ", ")
  )
}

r_packages <- extract_namespace_packages(r_files)
missing_r_packages <- setdiff(r_packages, r_declared)
if (length(missing_r_packages) > 0L) {
  record_failure(
    "missing-r-deps",
    paste(sort(missing_r_packages), collapse = ", ")
  )
}

test_packages <- extract_namespace_packages(test_files)
missing_test_packages <- setdiff(test_packages, c(r_declared, pkg_name))
if (length(missing_test_packages) > 0L) {
  record_failure(
    "missing-test-deps",
    paste(sort(missing_test_packages), collapse = ", ")
  )
}

if (length(failures) > 0L) {
  stop(
    paste(
      c("Package hardening static regression failed:", failures),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

message("Package hardening static regression passed.")
