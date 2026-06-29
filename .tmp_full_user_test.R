log_file <- "/Users/cxy/Desktop/lwdid_r/CRAN提交/user_full_test_log.txt"
sink(log_file, split = TRUE)
cat("====== lwdid EXHAUSTIVE USER TEST ======\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("R:", R.version.string, "\n\n")

# --- Install & Load ---
cat("== INSTALL ==\n")
install.packages("/Users/cxy/Desktop/lwdid_r/lwdid-r/lwdid_0.1.0.tar.gz", repos = NULL, type = "source")
library(lwdid)
cat("Version:", as.character(packageVersion("lwdid")), "\n\n")

# --- Load all datasets ---
cat("== DATASETS ==\n")
data("castle"); data("smoking"); data("walmart")
cat("castle:", nrow(castle), "x", ncol(castle), "\n")
cat("smoking:", nrow(smoking), "x", ncol(smoking), "\n")
cat("walmart:", nrow(walmart), "x", ncol(walmart), "\n")
cat("smoking cols:", paste(names(smoking), collapse = ", "), "\n\n")

# Helper
run_test <- function(desc, expr) {
  cat("--- ", desc, " ---\n")
  tryCatch({
    result <- eval(expr)
    if (inherits(result, "lwdid_result")) {
      cat("  ATT:", result$att, " SE:", result$se_att, " p:", result$pvalue, "\n")
      cat("  estimator:", result$estimator, " vce:", result$vce_type, "\n")
    } else {
      cat("  Result class:", paste(class(result), collapse = ", "), "\n")
    }
    cat("  STATUS: OK\n\n")
    invisible(result)
  }, warning = function(w) {
    cat("  WARNING:", conditionMessage(w), "\n")
    cat("  STATUS: OK (with warning)\n\n")
    suppressWarnings(eval(expr))
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    cat("  STATUS: FAILED\n\n")
    NULL
  })
}

# ======================================================================
# PART A: lwdid() - COMMON TIMING (smoking data)
# ======================================================================
cat("========================================\n")
cat("== PART A: lwdid() COMMON TIMING ==\n")
cat("========================================\n\n")

# A1: estimator options
cat("-- A1: estimator options --\n")
run_test("CT: estimator=ra", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ra")
))
run_test("CT: estimator=ipw (needs controls)", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ipw", controls="retprice")
))
run_test("CT: estimator=ipwra (needs controls)", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ipwra", controls="retprice")
))
run_test("CT: estimator=psm (needs controls)", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="psm", controls="retprice")
))

# A2: rolling options
cat("-- A2: rolling options --\n")
for (r in c("demean", "detrend")) {
  run_test(paste0("CT: rolling=", r), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling=.(r), estimator="ra")
  ))
}

# A3: aggregate options
cat("-- A3: aggregate options --\n")
for (a in c("overall", "cohort", "event_time")) {
  run_test(paste0("CT: aggregate=", a), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ra", aggregate=.(a))
  ))
}

# A4: vce options
cat("-- A4: vce options --\n")
for (v in c("homoskedastic", "HC0", "HC1", "HC2", "HC3", "HC4")) {
  run_test(paste0("CT: vce=", v), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ra", vce=.(v))
  ))
}
run_test("CT: vce=cluster", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ra", vce="cluster", cluster_var="state")
))

# A5: alpha options
cat("-- A5: alpha options --\n")
for (al in c(0.01, 0.05, 0.10)) {
  run_test(paste0("CT: alpha=", al), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ra", alpha=.(al))
  ))
}

# A6: ri options
cat("-- A6: ri options --\n")
run_test("CT: ri=TRUE, method=bootstrap", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ra", ri=TRUE, rireps=49L, seed=1L, ri_method="bootstrap")
))
run_test("CT: ri=TRUE, method=permutation", quote(
  lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ra", ri=TRUE, rireps=49L, seed=1L, ri_method="permutation")
))

# A7: WCB options
cat("-- A7: WCB options --\n")
for (wt in c("rademacher", "mammen", "webb")) {
  run_test(paste0("CT: WCB type=", wt), bquote(
    lwdid(data=smoking, y="cigsale", ivar="state", tvar="year", d="treat", post="post", rolling="demean", estimator="ra", vce="wild_cluster_bootstrap", cluster_var="state", wcb_reps=49L, wcb_type=.(wt), wcb_seed=1L)
  ))
}

# ======================================================================
# PART B: lwdid() - STAGGERED (castle data)
# ======================================================================
cat("========================================\n")
cat("== PART B: lwdid() STAGGERED ==\n")
cat("========================================\n\n")

# B1: control_group options
cat("-- B1: control_group options --\n")
for (cg in c("not_yet_treated", "never_treated")) {
  run_test(paste0("Stag: control_group=", cg), bquote(
    lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ra", aggregate="overall", control_group=.(cg))
  ))
}

# B2: estimators
cat("-- B2: estimator options --\n")
run_test("Stag: ra", quote(
  lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ra", aggregate="overall")
))
run_test("Stag: ipw (with controls, never_treated)", quote(
  lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ipw", aggregate="overall", control_group="never_treated", controls=c("police", "income"))
))
run_test("Stag: ipwra (with controls, never_treated)", quote(
  lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ipwra", aggregate="overall", control_group="never_treated", controls=c("police", "income"))
))

# B3: rolling
cat("-- B3: rolling --\n")
for (r in c("demean", "detrend")) {
  run_test(paste0("Stag: rolling=", r), bquote(
    lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling=.(r), estimator="ra", aggregate="overall")
  ))
}

# B4: aggregate
cat("-- B4: aggregate --\n")
for (a in c("overall", "cohort", "event_time")) {
  run_test(paste0("Stag: aggregate=", a), bquote(
    lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ra", aggregate=.(a))
  ))
}

# ======================================================================
# PART C: POST-ESTIMATION
# ======================================================================
cat("========================================\n")
cat("== PART C: POST-ESTIMATION ==\n")
cat("========================================\n\n")

res <- lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="demean", estimator="ra", aggregate="overall")

# C1: S3 methods
cat("-- C1: Standard S3 methods --\n")
run_test("coef()", quote(coef(res)))
run_test("confint()", quote(confint(res)))
run_test("confint(level=0.99)", quote(confint(res, level = 0.99)))
run_test("vcov()", quote(vcov(res)))
run_test("nobs()", quote(nobs(res)))
run_test("print()", quote(print(res)))
run_test("summary()", quote(summary(res)))

# C2: tidy/glance
cat("-- C2: tidy/glance --\n")
run_test("tidy(type=overall)", quote(tidy(res)))
run_test("tidy(type=overall, conf.int=FALSE)", quote(tidy(res, conf.int = FALSE)))
run_test("tidy(type=effects)", quote(tidy(res, type = "effects")))
run_test("tidy(type=event_time)", quote(tidy(res, type = "event_time")))
run_test("tidy(conf.level=0.99)", quote(tidy(res, conf.level = 0.99)))
run_test("glance()", quote(glance(res)))

# C3: fixest-style
cat("-- C3: fixest-style accessors --\n")
run_test("se()", quote(se(res)))
run_test("tstat()", quote(tstat(res)))
run_test("pvalue()", quote(pvalue(res)))
run_test("coeftable()", quote(coeftable(res)))

# C4: compare
cat("-- C4: compare() --\n")
res2 <- lwdid(data=castle, y="lhomicide", ivar="sid", tvar="year", gvar="gvar", rolling="detrend", estimator="ra", aggregate="overall")
run_test("compare(type=overall)", quote(compare(Demean=res, Detrend=res2, type="overall")))
run_test("compare(stars=FALSE)", quote(compare(Demean=res, Detrend=res2, stars=FALSE)))
run_test("compare(digits=5)", quote(compare(Demean=res, Detrend=res2, digits=5L)))

# C5: export
cat("-- C5: export functions --\n")
run_test("to_csv()", quote({f<-tempfile(fileext=".csv"); to_csv(res, file=f); cat("  wrote", file.size(f), "bytes\n"); f}))
run_test("to_dict()", quote(to_dict(res)))
run_test("to_latex()", quote(to_latex(res)))

# C6: plot
cat("-- C6: plot --\n")
run_test("plot()", quote({pdf(tempfile()); plot(res); dev.off(); "plotted"}))

cat("\n====== EXHAUSTIVE TEST COMPLETE ======\n")
cat("Total sections tested\n")
sink()
