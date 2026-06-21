files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in files) {
  tryCatch(
    parse(f),
    error = function(e) {
      cat("PARSE ERROR in", f, ":", conditionMessage(e), "\n")
    }
  )
}
cat("Done checking", length(files), "files\n")
