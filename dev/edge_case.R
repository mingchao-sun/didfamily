# One edge case per process, so that an R crash (segfault) is visible as a non-zero exit code.
#   Rscript dev/edge_case.R <case>     (jamovi-like library path, sandbox emulated)
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_estimators", "core_diag", "core_ddd", "core_dose", "core_panel", "core_code"))
    source(file.path(root, "R", paste0(f, ".R")))
source(file.path(root, "dev/sandbox.R")); dc_enable_sandbox()
suppressMessages(suppressWarnings(NULL))
case <- commandArgs(TRUE)[1]

make <- function(N = 300, T = 10, cohorts = c(4, 6, 0), years = NULL, seed = 1) {
    set.seed(seed)
    g <- sample(cohorts, N, TRUE)
    yrs <- if (is.null(years)) seq_len(T) else years
    d <- expand.grid(id = seq_len(N), t = seq_along(yrs)); d$year <- yrs[d$t]
    d$gi <- g[d$id]; d$g <- ifelse(d$gi > 0, yrs[pmin(d$gi, length(yrs))], 0)
    d$D <- as.integer(d$gi > 0 & d$t >= d$gi)
    d$y <- rnorm(N)[d$id] + 0.2 * d$t + d$D * 1 + rnorm(nrow(d), sd = .5)
    d[, c("id", "year", "y", "D", "g")]
}
run_all <- function(d, label) {
    out <- tryCatch({
        p <- dc_prepare(d, "y", "id", "year", "D", "indicator")
        rs <- list(dc_twfe(p, -3, 3), dc_cs(p, -3, 3), dc_sa(p, -3, 3), dc_bjs(p, -3, 3), dc_gardner(p, -3, 3), dc_stacked(p, -3, 3))
        sprintf("prepared; estimators ok: %s | messages: %s",
                paste(vapply(rs, function(r) if (isTRUE(r$ok)) "Y" else "n", ""), collapse = ""),
                paste(unique(vapply(rs, function(r) substr(r$message, 1, 70), "")[!vapply(rs, function(r) isTRUE(r$ok), TRUE)]), collapse = " / "))
    }, error = function(e) paste("prepare error:", substr(conditionMessage(e), 1, 110)))
    cat(sprintf("%-22s %s\n", label, out))
}

if (case == "no_treated")      { d <- make(cohorts = 0); d$D <- 0L; run_all(d, case) }
if (case == "all_same_cohort") { d <- make(cohorts = 4); run_all(d, case) }
if (case == "single_cohort")   { d <- make(cohorts = c(5, 0)); run_all(d, case) }
if (case == "tiny")            { d <- make(N = 12, T = 6, cohorts = c(3, 0)); run_all(d, case) }
if (case == "two_periods")     { d <- make(T = 2, cohorts = c(2, 0)); run_all(d, case) }
if (case == "constant_y")      { d <- make(); d$y <- 1; run_all(d, case) }
if (case == "treat_text")      { d <- make(); d$D <- ifelse(d$D == 1, "yes", "no"); run_all(d, case) }
if (case == "time_date")       { d <- make(); d$year <- sprintf("%d-01-01", 1999 + d$year); run_all(d, case) }
if (case == "time_from_zero")  { d <- make(years = 0:9); run_all(d, case) }
if (case == "time_negative")   { d <- make(years = -5:4); run_all(d, case) }
if (case == "time_step2")      { d <- make(years = seq(2000, 2018, by = 2)); run_all(d, case) }
if (case == "time_irregular")  { d <- make(T = 6, years = c(2000, 2001, 2003, 2007, 2008, 2012), cohorts = c(3, 5, 0)); run_all(d, case) }
if (case == "all_treated_late"){ d <- make(cohorts = c(10)); run_all(d, case) }
if (case == "huge_cohorts")    { d <- make(N = 400, T = 12, cohorts = c(2:11, 0)); run_all(d, case) }
if (case == "na_outcome")      { d <- make(); d$y[sample(nrow(d), 200)] <- NA; run_all(d, case) }
if (case == "inf_outcome")     { d <- make(); d$y[5] <- Inf; run_all(d, case) }
if (case == "one_unit")        { d <- make(N = 1, T = 8, cohorts = 4); run_all(d, case) }
if (case == "dup_rows")        { d <- make(); d <- rbind(d, d[1:5, ]); run_all(d, case) }
if (case == "factor_id")       { d <- make(); d$id <- paste0("unit_", d$id); run_all(d, case) }
if (case == "diag_edge") {     # diagnostics + DDD + 2x2 + panel on the same odd data
    d <- make(cohorts = c(5, 0)); p <- dc_prepare(d, "y", "id", "year", "D", "indicator")
    r1 <- dc_bacon(p); r2 <- dc_twfe_weights(p); r3 <- dc_pretests(p, -3); r4 <- tryCatch(dc_honest_prepare(p, "sa", -3, 3), error = function(e) conditionMessage(e))
    cat(sprintf("%-22s bacon ok=%s weights ok=%s pretests rows=%d honest=%s\n", case, r1$ok, r2$ok, nrow(r3), if (is.list(r4)) "ok" else substr(r4, 1, 60)))
}
if (case == "panel_all") {
    for (nm in c("no_treated", "all_same_cohort", "tiny", "time_date", "dup_rows")) {
        d <- switch(nm, no_treated = {x <- make(cohorts = 0); x$D <- 0L; x}, all_same_cohort = make(cohorts = 4), tiny = make(N = 12, T = 6, cohorts = c(3, 0)),
                    time_date = {x <- make(); x$year <- sprintf("%d-01-01", 1999 + x$year); x}, dup_rows = {x <- make(); rbind(x, x[1:5, ])})
        r <- tryCatch(dc_panel_check(d, "id", "year", "D", "indicator", dep = "y"), error = function(e) conditionMessage(e))
        cat(sprintf("panel_check %-16s %s\n", nm, if (is.list(r)) paste(sum(r$checks$status %in% c("problem", "warn")), "problem/warn rows") else paste("ERROR", substr(r, 1, 80))))
    }
}
