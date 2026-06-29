log_file <- "/Users/cxy/Desktop/lwdid_r/CRAN提交/user_exhaustive_test_log.txt"
con <- file(log_file, "w")
sink(con, type = "output")
sink(con, type = "message")

cat("====== lwdid EXHAUSTIVE PARAMETER TEST ======\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("R:", R.version.string, "\n\n")

library(lwdid)
data("castle"); data("smoking"); data("walmart")

run <- function(desc, expr) {
  cat(sprintf("[TEST] %s\n", desc))
  t0 <- proc.time()
  ok <- TRUE
  result <- tryCatch(
    withCallingHandlers(
      eval(expr, envir = parent.frame()),
      warning = function(w) {
        cat(sprintf("  WARNING: %s\n", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      elapsed <- (proc.time() - t0)[3]
      cat(sprintf("  ERROR: %s (%.1fs)\n", conditionMessage(e), elapsed))
      cat("  STATUS: FAILED\n\n")
      ok <<- FALSE
      NULL
    }
  )
  if (ok) {
    elapsed <- (proc.time() - t0)[3]
    r <- result
    if (inherits(r, "lwdid_result")) {
      cat(sprintf("  -> ATT=%.6f SE=%.6f p=%.4f vce=%s (%.1fs)\n", r$att, r$se_att, r$pvalue, r$vce_type, elapsed))
    } else if (is.data.frame(r)) {
      cat(sprintf("  -> data.frame %dx%d (%.1fs)\n", nrow(r), ncol(r), elapsed))
    } else {
      cat(sprintf("  -> %s (%.1fs)\n", paste(class(r), collapse=","), elapsed))
    }
    cat("  STATUS: OK\n\n")
  }
  invisible(result)
}

# ========================================================
cat("\n######## SECTION 1: COMMON TIMING — ALL ESTIMATORS ########\n\n")

for (est in c("ra", "ipw", "ipwra", "psm")) {
  if (est == "ra") {
    run(paste0("CT estimator=", est), bquote(
      lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator=.(est))
    ))
  } else {
    run(paste0("CT estimator=", est, " (with controls)"), bquote(
      lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator=.(est), controls="retprice")
    ))
  }
}

# ========================================================
cat("\n######## SECTION 2: CT — ALL ROLLING ########\n\n")

for (r in c("demean", "detrend")) {
  run(paste0("CT rolling=", r), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling=.(r), estimator="ra")
  ))
}

# ========================================================
cat("\n######## SECTION 3: CT — ALL AGGREGATE ########\n\n")

for (a in c("overall", "cohort", "event_time")) {
  run(paste0("CT aggregate=", a), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator="ra", aggregate=.(a))
  ))
}

# ========================================================
cat("\n######## SECTION 4: CT — ALL VCE OPTIONS ########\n\n")

for (v in c("hc0", "hc1", "hc2", "hc3", "hc4", "HC0", "HC1", "HC2", "HC3", "HC4",
            "robust", "Robust", "ROBUST",
            "homoskedastic", "Homoskedastic",
            "cluster")) {
  if (v %in% c("cluster")) {
    run(paste0("CT vce=", v), bquote(
      lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator="ra", vce=.(v), cluster_var="state")
    ))
  } else {
    run(paste0("CT vce=", v), bquote(
      lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator="ra", vce=.(v))
    ))
  }
}

# WCB aliases
for (v in c("bootstrap", "wcb", "wild_cluster_bootstrap")) {
  run(paste0("CT vce=", v), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator="ra", vce=.(v), cluster_var="state", wcb_reps=29L, wcb_seed=1L)
  ))
}

# ========================================================
cat("\n######## SECTION 5: CT — ALPHA ########\n\n")

for (al in c(0.01, 0.05, 0.10, 0.20)) {
  run(paste0("CT alpha=", al), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator="ra", alpha=.(al))
  ))
}

# ========================================================
cat("\n######## SECTION 6: CT — RI ########\n\n")

for (m in c("bootstrap", "permutation")) {
  run(paste0("CT ri method=", m), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator="ra", ri=TRUE, rireps=29L, seed=1L, ri_method=.(m))
  ))
}

# ========================================================
cat("\n######## SECTION 7: CT — WCB WEIGHT TYPES ########\n\n")

for (wt in c("rademacher", "mammen", "webb")) {
  run(paste0("CT wcb_type=", wt), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="demean", estimator="ra", vce="bootstrap", cluster_var="state", wcb_reps=29L, wcb_type=.(wt), wcb_seed=1L)
  ))
}

# ========================================================
cat("\n######## SECTION 8: STAGGERED — CONTROL GROUPS ########\n\n")

for (cg in c("not_yet_treated", "never_treated")) {
  run(paste0("Stag control_group=", cg), bquote(
    lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ra", aggregate="overall", control_group=.(cg))
  ))
}

# ========================================================
cat("\n######## SECTION 9: STAGGERED — ALL ESTIMATORS ########\n\n")

run("Stag ra", quote(
  lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ra", aggregate="overall", control_group="never_treated")
))
run("Stag ipw", quote(
  lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ipw", aggregate="overall", control_group="never_treated", controls=c("police","income"))
))
run("Stag ipwra", quote(
  lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ipwra", aggregate="overall", control_group="never_treated", controls=c("police","income"))
))

# ========================================================
cat("\n######## SECTION 10: STAGGERED — ROLLING ########\n\n")

for (r in c("demean", "detrend")) {
  run(paste0("Stag rolling=", r), bquote(
    lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling=.(r), estimator="ra", aggregate="overall")
  ))
}

# ========================================================
cat("\n######## SECTION 11: STAGGERED — AGGREGATE ########\n\n")

for (a in c("overall", "cohort", "event_time")) {
  run(paste0("Stag aggregate=", a), bquote(
    lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ra", aggregate=.(a))
  ))
}

# ========================================================
cat("\n######## SECTION 12: POST-ESTIMATION ########\n\n")

res <- lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ra", aggregate="overall")

run("print()", quote(print(res)))
run("summary()", quote(summary(res)))
run("coef()", quote(coef(res)))
run("confint()", quote(confint(res)))
run("confint(level=0.99)", quote(confint(res, level=0.99)))
run("confint(level=0.90)", quote(confint(res, level=0.90)))
run("vcov()", quote(vcov(res)))
run("nobs()", quote(nobs(res)))

# ========================================================
cat("\n######## SECTION 13: TIDY/GLANCE ########\n\n")

run("tidy()", quote(tidy(res)))
run("tidy(type='overall')", quote(tidy(res, type="overall")))
run("tidy(type='effects')", quote(tidy(res, type="effects")))
run("tidy(type='event_time')", quote(tidy(res, type="event_time")))
run("tidy(conf.int=FALSE)", quote(tidy(res, conf.int=FALSE)))
run("tidy(conf.level=0.99)", quote(tidy(res, conf.level=0.99)))
run("tidy(conf.level=0.90)", quote(tidy(res, conf.level=0.90)))
run("glance()", quote(glance(res)))

# ========================================================
cat("\n######## SECTION 14: FIXEST-STYLE ########\n\n")

run("se()", quote(se(res)))
run("tstat()", quote(tstat(res)))
run("pvalue()", quote(pvalue(res)))
run("coeftable()", quote(coeftable(res)))

# ========================================================
cat("\n######## SECTION 15: COMPARE ########\n\n")

res2 <- lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="detrend", estimator="ra", aggregate="overall")

run("compare(type='overall')", quote(compare(Demean=res, Detrend=res2, type="overall")))
run("compare(type='effects')", quote(compare(Demean=res, Detrend=res2, type="effects")))
run("compare(stars=FALSE)", quote(compare(Demean=res, Detrend=res2, stars=FALSE)))
run("compare(stars=TRUE)", quote(compare(Demean=res, Detrend=res2, stars=TRUE)))
run("compare(digits=2)", quote(compare(Demean=res, Detrend=res2, digits=2L)))
run("compare(digits=5)", quote(compare(Demean=res, Detrend=res2, digits=5L)))

# ========================================================
cat("\n######## SECTION 16: EXPORT ########\n\n")

run("to_csv()", quote({f<-tempfile(fileext=".csv"); to_csv(res, file=f); cat(sprintf("  wrote %d bytes\n", file.size(f))); f}))
run("to_dict()", quote({d<-to_dict(res); cat(sprintf("  list length=%d\n", length(d))); d}))
run("to_latex()", quote(to_latex(res)))

# ========================================================
cat("\n######## SECTION 17: PLOT ########\n\n")

run("plot(res)", quote({pdf(tf<-tempfile(fileext=".pdf")); plot(res); dev.off(); cat(sprintf("  %d bytes\n", file.size(tf))); tf}))

# ========================================================
cat("\n######## SECTION 18: ERROR SCENARIOS ########\n\n")

run("ERROR: missing y column", quote(
  lwdid(data=smoking, y="nonexist", ivar="state", tvar="year", d="d", post="post")
))
run("ERROR: invalid estimator", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", estimator="invalid")
))
run("ERROR: invalid rolling", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", rolling="invalid")
))
run("ERROR: ipw without controls", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", estimator="ipw")
))
run("ERROR: invalid vce", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", vce="nonsense")
))
run("ERROR: cluster without cluster_var", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="d", post="post", vce="cluster")
))

# ========================================================
cat("\n######## SECTION 19: WALMART DATA ########\n\n")

run("walmart staggered demean ra", quote(
  lwdid(data=walmart, y="log_retail_emp", ivar="fips", tvar="year", gvar="g", rolling="demean", estimator="ra", aggregate="overall")
))
run("walmart staggered detrend ra", quote(
  lwdid(data=walmart, y="log_retail_emp", ivar="fips", tvar="year", gvar="g", rolling="detrend", estimator="ra", aggregate="overall")
))

# ========================================================
cat("\n######## SECTION 20: CITATION ########\n\n")

run("citation()", quote(citation("lwdid")))

cat("\n====== ALL TESTS COMPLETE ======\n")
sink(type = "message")
sink(type = "output")
close(con)
