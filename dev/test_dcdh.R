# The dynamic DID estimators for non-absorbing treatments (R/core_dcdh.R) against the reference implementation
# DIDmultiplegtDYN (1.x/2.1 releases need no polars). Set DCDH_REFLIB to a library that holds it, else only
# internal checks run.   Run from the module root:  Rscript dev/test_dcdh.R
root <- if (file.exists("R/core_dcdh.R")) "." else ".."
source(file.path(root, "R", "core_dcdh.R"))
lib <- Sys.getenv("DCDH_REFLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
have_ref <- requireNamespace("DIDmultiplegtDYN", quietly = TRUE)
if (have_ref) suppressMessages(library(DIDmultiplegtDYN)) else cat("DIDmultiplegtDYN not available: only internal checks\n")
nfail <- 0
chk <- function(name, got, want, tol) {
    d <- max(abs(got - want), na.rm = TRUE); ok <- is.finite(d) && d < tol
    cat(sprintf("%-58s %s  max |diff| = %.1e\n", name, if (ok) "ok  " else "FAIL", d)); if (!ok) nfail <<- nfail + 1
}
make <- function(seed, N = 200, Tn = 10, p_treat = 0.6, p_base1 = 0, p_off = 0.4, drop = 0, never = TRUE) {
    set.seed(seed)
    d <- expand.grid(id = 1:N, t = 1:Tn)
    base <- rbinom(N, 1, p_base1)
    first <- sample(3:8, N, TRUE)
    sw <- if (never) ifelse(runif(N) < p_treat, first, NA) else first
    back <- ifelse(!is.na(sw) & runif(N) < p_off, sw + sample(1:2, N, TRUE), NA)      # return to baseline
    on <- !is.na(sw[d$id]) & d$t >= sw[d$id] & (is.na(back[d$id]) | d$t < back[d$id])
    d$D <- ifelse(on, 1 - base[d$id], base[d$id])
    d$y <- rnorm(N)[d$id] + 0.1 * d$t + 1.5 * ifelse(d$D == 1, 1, 0) + rnorm(nrow(d))
    if (drop > 0) d <- d[runif(nrow(d)) > drop, ]
    d
}
cases <- list(
    "switching on and off, never-treated present" = list(seed = 1),
    "baseline 0 and 1 (switching in and out)"     = list(seed = 2, p_base1 = 0.4),
    "unbalanced panel (6% rows missing)"          = list(seed = 3, drop = 0.06),
    "everyone eventually switches (no never-treated)" = list(seed = 4, never = FALSE, p_off = 0.2),
    "absorbing treatment only"                    = list(seed = 5, p_off = 0))
for (nm in names(cases)) {
    cat("\n---", nm, "\n")
    d <- do.call(make, cases[[nm]])
    me <- dc_dcdh(d$y, d$id, d$t, d$D, effects = 3, placebo = 2)
    chk("all estimates finite", as.numeric(is.finite(me$table$estimate)) , 1, 0.5)
    if (have_ref) {
        ref <- suppressWarnings(did_multiplegt_dyn(df = d, outcome = "y", group = "id", time = "t", treatment = "D", effects = 3, placebo = 2, graph_off = TRUE))
        rt <- rbind(ref$results$Effects, ref$results$Placebos)
        chk("estimates", me$table$estimate, as.numeric(ref$coef$b), 1e-5)     # the reference rounds its estimates to 5 decimals
        chk("standard errors", me$table$se, sqrt(diag(ref$coef$vcov)), 1e-8)
        chk("number of switchers", me$table$n_switchers, rt[, 6], 0.5)
        jp <- dc_dcdh_joint(me, "placebo")
        chk("joint placebo p-value", jp$p, ref$results$p_jointplacebo, 1e-5)
    }
}
cat(if (nfail == 0) "\nAll dynamic DID checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
