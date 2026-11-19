# The module's dosedid and dddcompare against the R packages they wrap (contdid, triplediff) called directly on the raw data.
# Checks that the module passes the data and options through unchanged. Needs both packages (library named by DCDH_REFLIB). Run from didfamily/.
lib <- Sys.getenv("DCDH_REFLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
if (!requireNamespace("contdid", quietly = TRUE) || !requireNamespace("triplediff", quietly = TRUE)) { cat("contdid or triplediff not installed: skipped\n"); quit(status = 0) }
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_dose", "core_ddd", "core_estimators")) source(file.path(root, "R", paste0(f, ".R")))
nfail <- 0
cmp <- function(label, a, b, tol = 1e-8) { ok <- isTRUE(all.equal(unname(a), unname(b), tolerance = tol)); cat(sprintf("[%s] %-46s module %s  package %s\n", if (ok) "PASS" else "FAIL", label, format(signif(a[1], 8)), format(signif(b[1], 8)))); if (!ok) nfail <<- nfail + 1 }

# ---- continuous treatment: contdid ----
# contdid draws a multiplier bootstrap even with bstrap = FALSE (it sets the overall s.e.), the module fixes seed 1234 and so does this test before the direct call.
set.seed(4); N <- 800; T <- 4
g <- sample(c(3, 0), N, TRUE, prob = c(.6, .4)); dose <- ifelse(g > 0, runif(N, .1, 1), 0)
dat <- expand.grid(id = 1:N, t = 1:T); dat$g <- g[dat$id]; dat$dose <- dose[dat$id]
dat$y <- rnorm(N)[dat$id] + .3 * dat$t + ifelse(dat$g > 0 & dat$t >= dat$g, 2 * dat$dose + dat$dose^2, 0) + rnorm(nrow(dat), sd = .3)
d <- dc_prepare(dat, "y", "id", "t", "g", "cohort", dose = "dose")
for (cfg in list(list(deg = 1, kn = 0, ctl = "nevertreated"), list(deg = 2, kn = 0, ctl = "nevertreated"), list(deg = 3, kn = 1, ctl = "nevertreated"), list(deg = 2, kn = 0, ctl = "notyettreated"))) {
    m <- dc_dose(d, cfg$ctl, cfg$deg, cfg$kn); m2 <- dc_dose(d, cfg$ctl, cfg$deg, cfg$kn)
    cmp("dose result is reproducible", m$att_se, m2$att_se, 0)
    set.seed(1234)
    p <- contdid::cont_did(yname = "y", dname = "dose", gname = "g", tname = "t", idname = "id", data = dat, target_parameter = "level", aggregation = "dose",
                           treatment_type = "continuous", dose_est_method = "parametric", degree = cfg$deg, num_knots = cfg$kn, control_group = cfg$ctl, bstrap = FALSE, cband = FALSE, print_details = FALSE)
    lab <- sprintf("dose degree %d knots %d, %s", cfg$deg, cfg$kn, cfg$ctl)
    cmp(paste(lab, "overall ATT"), m$att, p$overall_att); cmp(paste(lab, "overall ATT s.e."), m$att_se, p$overall_att_se)
    cmp(paste(lab, "overall ACRT"), m$acrt, p$overall_acrt); cmp(paste(lab, "dose-response curve"), m$curve$att, as.numeric(p$att.d))
}

# ---- triple differences: triplediff ----
dd <- read.csv(file.path(root, "dev/csv/sim_ddd.csv")); dp <- dc_prepare(dd, "y", "id", "year", "g", "cohort", part = "elig")
for (cfg in list(list("dr", "nevertreated"), list("reg", "nevertreated"), list("ipw", "nevertreated"))) {
    m <- dc_ddd(dp, -3, 3, cfg[[2]], cfg[[1]])
    fit <- triplediff::ddd(yname = "y", tname = "year", idname = "id", gname = "g", pname = "elig", xformla = ~1, data = dd, control_group = cfg[[2]], base_period = "universal", est_method = cfg[[1]])
    s <- triplediff::agg_ddd(fit, type = "simple")$aggte_ddd
    lab <- sprintf("DDD %s", cfg[[1]])
    cmp(paste(lab, "overall ATT"), m$att, s$overall.att); cmp(paste(lab, "overall ATT s.e."), m$se, s$overall.se)
}
cat(if (nfail == 0) "\nAll wrapper checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
