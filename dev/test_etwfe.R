# The native extended TWFE (R/core_etwfe.R) against the R package etwfe: coefficients and the four aggregations (simple, group, calendar, event) with their
# standard errors, linear, Poisson and logit, with and without covariates, both control groups. Needs etwfe (default library or DCDH_REFLIB).
# Run from the module root:  Rscript dev/test_etwfe.R
root <- if (file.exists("R/core_etwfe.R")) "." else ".."
lib <- Sys.getenv("DCDH_REFLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(jmvcore))
for (f in c("core_prepare", "core_estimators", "core_etwfe")) source(file.path(root, "R", paste0(f, ".R")))
if (!requireNamespace("etwfe", quietly = TRUE)) { cat("etwfe is not installed: skipped\n"); quit(status = 0) }
nfail <- 0
cmp <- function(label, a, b, tol) { ok <- length(a) == length(b) && all(is.finite(a)) && max(abs(a - b)) < tol; cat(sprintf("[%s] %-70s max diff %.1e\n", if (ok) "PASS" else "FAIL", label, if (length(a) == length(b)) max(abs(a - b)) else NA)); if (!ok) nfail <<- nfail + 1 }
ok <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) nfail <<- nfail + 1 }

mp <- read.csv(file.path(root, "dev/csv/mpdta.csv"))
sim <- function(seed, N = 400, Tn = 8, family = "gaussian") {
    set.seed(seed)
    d <- expand.grid(id = 1:N, year = 1:Tn); first <- sample(c(4, 5, 6, 0), N, TRUE, c(.25, .25, .2, .3))
    x1 <- rnorm(N); x2 <- rbinom(N, 1, .4); d$g <- first[d$id]; d$x1 <- x1[d$id]; d$x2 <- x2[d$id]
    treated <- d$g > 0 & d$year >= d$g
    eff <- (0.4 + 0.1 * (d$year - d$g)) * treated * (1 + 0.5 * d$x1)
    lin <- 0.3 * d$x1 + 0.2 * d$x2 + 0.1 * d$year + rnorm(N)[d$id] * 0.3
    if (family == "gaussian") d$y <- lin + eff + rnorm(nrow(d))
    else if (family == "poisson") d$y <- rpois(nrow(d), exp(0.5 + 0.2 * lin + 0.5 * eff))
    else d$y <- rbinom(nrow(d), 1, plogis(0.2 * lin + 0.8 * eff - 0.3))
    d
}
prep <- function(dat, y, id, tm, g, covs = NULL) dc_prepare(dat, y, id, tm, g, "cohort", covs = covs)

check <- function(label, dat, ycol, idcol, tcol, gcol, control, family, covs = NULL, tol = 2e-6) {
    d <- prep(dat, ycol, idcol, tcol, gcol, covs)
    cx <- d$info$covs; cc <- if (is.null(cx) || !nrow(cx)) character(0) else cx$code
    mine <- dc_etwfe_fit(d$data, cc, control, family)
    df <- d$data
    fml <- if (length(cc)) stats::as.formula(paste("dc_y ~", paste(cc, collapse = " + "))) else dc_y ~ 0
    gr <- if (any(df$dc_g == 0)) 0 else NULL
    ref <- suppressWarnings(suppressMessages(etwfe::etwfe(fml, tvar = dc_t, gvar = dc_g, data = df, gref = gr, cgroup = if (control == "nevertreated") "never" else "notyet",
                                                          family = if (family == "gaussian") NULL else if (family == "poisson") "poisson" else "logit", vcov = ~dc_cl)))
    for (ty in c("simple", "group", "calendar", "event")) {
        a <- dc_etwfe_agg(mine, ty)
        r <- suppressWarnings(suppressMessages(etwfe::emfx(ref, type = ty, post_only = (ty != "event" || control == "notyettreated"))))
        r <- as.data.frame(r)
        key <- switch(ty, simple = rep(1, nrow(r)), group = r$dc_g, calendar = r$dc_t, event = r$event)
        r <- r[order(key), ]; keyr <- sort(key)
        if (ty == "event") {                               # etwfe also lists reference periods whose effect is exactly zero
            extra <- !(keyr %in% a$key); ok(all(abs(r$estimate[extra]) < 1e-12), sprintf("%s: event times only in etwfe are zero effects (%s)", label, paste(keyr[extra], collapse = ",")))
            r <- r[!extra, ]; keyr <- keyr[!extra]
        }
        ok_keys <- identical(as.numeric(a$key), as.numeric(keyr))
        if (!ok_keys) { ok(FALSE, sprintf("%s, %s: same %s groups (mine %s / etwfe %s)", label, ty, ty, paste(a$key, collapse = ","), paste(keyr, collapse = ","))); next }
        cmp(sprintf("%s: %s effects and standard errors", label, ty), c(a$est, a$se), c(r$estimate, r$std.error), tol)
    }
}
check("mpdta, linear, not yet treated", mp, "lemp", "countyreal", "year", "first.treat", "notyettreated", "gaussian")
check("mpdta, linear, never treated", mp, "lemp", "countyreal", "year", "first.treat", "nevertreated", "gaussian")
check("mpdta, linear, covariate lpop", mp, "lemp", "countyreal", "year", "first.treat", "notyettreated", "gaussian", "lpop")
check("mpdta, linear, covariate lpop, never treated", mp, "lemp", "countyreal", "year", "first.treat", "nevertreated", "gaussian", "lpop")
sg <- sim(1, family = "gaussian")
check("simulated, linear, two covariates", sg, "y", "id", "year", "g", "notyettreated", "gaussian", c("x1", "x2"))
sp <- sim(2, family = "poisson")
check("simulated, Poisson, not yet treated", sp, "y", "id", "year", "g", "notyettreated", "poisson")
check("simulated, Poisson, covariates, never treated", sp, "y", "id", "year", "g", "nevertreated", "poisson", c("x1", "x2"), tol = 5e-6)
sl <- sim(3, family = "logit")
check("simulated, logit, not yet treated", sl, "y", "id", "year", "g", "notyettreated", "logit")
check("simulated, logit, covariate, never treated", sl, "y", "id", "year", "g", "nevertreated", "logit", "x1", tol = 5e-6)
# without never-treated units: the last cohort is the control
nn <- sg[sg$g != 0, ]; nn$g[nn$g == 6] <- 6
check("simulated, linear, no never-treated units", nn, "y", "id", "year", "g", "notyettreated", "gaussian")
# the adapter and the guards
d <- prep(mp, "lemp", "countyreal", "year", "first.treat")
a <- dc_etwfe(d, -4, 4, "nevertreated", "gaussian")
ok(a$ok && is.finite(a$att) && any(a$dynamic$e < -1), sprintf("adapter: overall ATT %.4f (s.e. %.4f) and an event study with leads", a$att, a$se))
b <- dc_bjs(d, -4, 4)
ok(abs(a$att - b$att) < 0.1, sprintf("without covariates the overall ETWFE effect %.4f is close to the imputation estimate %.4f", dc_etwfe(d, -4, 4, "notyettreated")$att, b$att))
ok(!dc_etwfe(prep(nn, "y", "id", "year", "g"), -4, 4, "nevertreated")$ok, "never-treated control without never-treated units is refused")
ok(!dc_etwfe(prep(sg, "y", "id", "year", "g"), -4, 4, "notyettreated", "logit")$ok, "the logit model needs a 0/1 outcome")
ok(!dc_etwfe(prep(sg, "y", "id", "year", "g"), -4, 4, "notyettreated", "poisson")$ok || all(sg$y >= 0), "the Poisson model needs a non-negative outcome")
cat(if (nfail == 0) "\nAll extended TWFE checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
