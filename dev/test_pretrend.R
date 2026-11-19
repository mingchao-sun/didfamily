# Power of the pre-trend test (R/core_pretrend.R): closed forms, Monte Carlo, and an end-to-end check against the rejection rate of the test in simulated panels.
# Run from the module root:  Rscript dev/test_pretrend.R
root <- if (file.exists("R/core_pretrend.R")) "." else ".."
lib <- Sys.getenv("DCDH_REFLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(jmvcore))
for (f in c("core_prepare", "core_estimators", "core_diag", "core_pretrend")) source(file.path(root, "R", paste0(f, ".R")))
nfail <- 0
ok <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) nfail <<- nfail + 1 }
crit <- qnorm(0.975)

# ---- closed forms
se <- 0.2
for (d in c(0, 0.2, 0.5)) ok(abs(dc_pt_reject(d, se^2, crit) - (1 - (pnorm(crit - d / se) - pnorm(-crit - d / se)))) < 1e-12, sprintf("one pre-period, shift %.1f: power equals the normal formula", d))
S <- diag(c(0.1, 0.15, 0.2)^2); dl <- c(0.1, 0.25, 0.3)
ok(abs(dc_pt_reject(dl, S, crit) - (1 - prod(pnorm(crit - dl / sqrt(diag(S))) - pnorm(-crit - dl / sqrt(diag(S)))))) < 1e-6, "independent coefficients: power equals the product formula")
ok(abs(dc_pt_reject(rep(0, 4), diag(4) * 0.1^2, crit) - (1 - 0.95^4)) < 1e-6, "size with four independent pre-periods is 1 - 0.95^4 = 18.5%")
# correlated coefficients against Monte Carlo
set.seed(1); A <- matrix(rnorm(9), 3); Sg <- crossprod(A) * 0.01 + diag(3) * 0.01; dl2 <- c(-0.2, 0.15, 0.3)
z <- MASS::mvrnorm(2e6, dl2, Sg); rej <- mean(apply(abs(z) > rep(crit * sqrt(diag(Sg)), each = nrow(z)), 1, any))
ok(abs(dc_pt_reject(dl2, Sg, crit) - rej) < 2e-3, sprintf("correlated coefficients: power %.4f against Monte Carlo %.4f", dc_pt_reject(dl2, Sg, crit), rej))
# slope for a target power
for (tg in c(0.5, 0.8)) { sl <- dc_pt_slope(S, c(-4, -3, -2), -1, tg, crit); ok(abs(dc_pt_reject((c(-4, -3, -2) + 1) * sl, S, crit) - tg) < 1e-5, sprintf("slope for %.0f%% power reproduces that power (slope %.4f)", 100 * tg, sl)) }
ok(dc_pt_slope(S, c(-4, -3, -2), -1, 0.8, crit) > dc_pt_slope(S, c(-4, -3, -2), -1, 0.5, crit), "a larger power needs a larger slope")
ok(abs(dc_pt_reject(-dl, S, crit) - dc_pt_reject(dl, S, crit)) < 1e-6, "power is symmetric in the sign of the violation")

# ---- end to end on a simulated staggered panel with a differential linear pre-trend
sim <- function(seed, slope, N = 400, Tn = 9) {
    set.seed(seed)
    d <- expand.grid(id = 1:N, t = 1:Tn); first <- sample(c(5, 7, 0), N, TRUE, c(.35, .3, .35)); g <- first[d$id]
    ever <- g > 0
    d$y <- rnorm(N)[d$id] + 0.2 * d$t + slope * ever[d$id] * (d$t - 4) + 1.0 * (ever[d$id] & d$t >= g) + rnorm(nrow(d), 0, 1)
    d$g <- g
    d
}
prep_of <- function(dat) { dd <- dc_prepare(dat, "y", "id", "t", "g", "cohort"); dc_honest_prepare(dd, "sa", -4, 4) }
p0 <- prep_of(sim(1, 0.15)); res <- dc_pretrend(p0, slope = 0.15)
ok(length(res$e_pre) >= 2 && is.finite(res$power) && res$power > res$size, sprintf("end to end: %d pre-periods, power against the trend 0.15 is %.2f, size %.2f", length(res$e_pre), res$power, res$size))
# the predicted power agrees with the rejection rate of the test over repeated samples
rejs <- replicate(80, { pr <- prep_of(sim(sample.int(1e6, 1), 0.15)); any(abs(pr$beta[pr$e < -1] / sqrt(diag(pr$sigma)[pr$e < -1])) > crit) })
ok(abs(mean(rejs) - res$power) < 0.15, sprintf("predicted power %.2f against the rejection rate %.2f over 80 simulated panels", res$power, mean(rejs)))
rej0 <- replicate(80, { pr <- prep_of(sim(sample.int(1e6, 1), 0)); any(abs(pr$beta[pr$e < -1] / sqrt(diag(pr$sigma)[pr$e < -1])) > crit) })
ok(abs(mean(rej0) - res$size) < 0.12, sprintf("with parallel trends the test rejects in %.2f of panels, predicted size %.2f", mean(rej0), res$size))
r1 <- dc_pretrend(p0, slope = NA); ok(is.na(r1$power) && all(r1$slope_for_power > 0) && r1$slope_for_power[2] > r1$slope_for_power[1], "without a hypothesized slope only the slopes for 50% and 80% power are returned")
ok(inherits(try(dc_pretrend(list(beta = 1, sigma = matrix(1), e = 0)), silent = TRUE), "try-error"), "no pre-treatment periods is refused")
cat(if (nfail == 0) "\nAll pre-trend power checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
