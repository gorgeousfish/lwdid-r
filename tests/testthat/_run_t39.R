devtools::load_all(".")

make_panel <- function(N = 8L, TT = 6L, S = 4L, n_treated = 5L,
                       tau = 3.0, sd_noise = 0.5, seed = 123L,
                       add_controls = FALSE) {
  alpha <- seq(10, by = 10, length.out = N)
  set.seed(seed)
  dt <- data.table::CJ(id = seq_len(N), time = seq_len(TT))
  dt[, d := as.integer(id <= n_treated)]
  dt[, post := as.integer(time >= S)]
  dt[, y := alpha[id] + tau * d * post + rnorm(.N, sd = sd_noise)]
  dt
}

quiet_lwdid <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

tryCatch({
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          vce = "hc1")
  )
  cat("class:", paste(class(result), collapse = ", "), "\n")
  cat("att:", result$att, "\n")
  cat("se:", result$se_att, "\n")
  cat("vce_type:", result$vce_type, "\n")
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("class:", paste(class(e), collapse = ", "), "\n")
}, warning = function(w) {
  cat("WARNING:", conditionMessage(w), "\n")
  cat("class:", paste(class(w), collapse = ", "), "\n")
})
