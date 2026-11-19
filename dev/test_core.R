# Core-layer test. Run in a jamovi-like library environment:
#   R_LIBS_USER=notthere R_LIBS="<module>/R;<jamovi>/Resources/modules/jmv/R;<jamovi>/Resources/modules/base/R"
#   Rscript dev/test_core.R
root <- if (file.exists("R/core_prepare.R")) "." else ".."
source(file.path(root, "R/core_prepare.R"))
source(file.path(root, "R/core_estimators.R")); source(file.path(root, "R/core_did2s.R")); source(file.path(root, "R/vendor_honestdid.R")); source(file.path(root, "R/core_diag.R")); source(file.path(root, "R/core_ddd.R"))

check <- function(cond, msg) {
    cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
    if (!isTRUE(cond)) quit(status = 1)
}

## ---- simulated staggered panel with heterogeneous, dynamic effects ----
set.seed(1)
N <- 600; T <- 12
cohort <- sample(c(5, 7, 9, 0), N, replace = TRUE, prob = c(.25, .25, .25, .25))
slope  <- c(`5` = 0.2, `7` = 0.6, `9` = 1.0)    # effect grows with exposure, by cohort
sim <- expand.grid(id = 1:N, year = 1:T)
sim$g <- cohort[sim$id]
sim$tau <- ifelse(sim$g > 0 & sim$year >= sim$g,
                  slope[as.character(sim$g)] * (sim$year - sim$g + 1), 0)
sim$tau[is.na(sim$tau)] <- 0
alpha <- rnorm(N); lambda <- rnorm(T, sd = .5)
sim$y <- alpha[sim$id] + lambda[sim$year] + sim$tau + rnorm(nrow(sim), sd = .5)
sim$D <- as.integer(sim$g > 0 & sim$year >= sim$g)
true_att <- mean(sim$tau[sim$D == 1])
cat(sprintf("true mean effect on treated unit-periods: %.3f\n", true_att))

d <- dc_prepare(sim, "y", "id", "year", "D", "indicator")
check(d$info$n_units == N && d$info$balanced, "prepare: units and balance")
check(identical(sort(as.numeric(names(d$info$cohorts))), c(5, 7, 9)), "prepare: cohorts recovered from indicator")

d2 <- dc_prepare(sim, "y", "id", "year", "g", "cohort")
check(all(d$data$dc_d == d2$data$dc_d), "prepare: cohort input equals indicator input")

tw <- dc_twfe(d, -4, 4); cs <- dc_cs(d, -4, 4, "nevertreated"); sa <- dc_sa(d, -4, 4)
cs2 <- dc_cs(d, -4, 4, "notyettreated")
for (r in list(tw, cs, sa, cs2))
    cat(sprintf("%-20s ATT %.3f (se %.3f) ok=%s %s\n", r$method, r$att, r$se, r$ok, r$message))

check(all(vapply(list(tw, cs, sa, cs2), function(r) r$ok, TRUE)), "all estimators ran")
check(abs(cs$att - true_att) < 0.15, "CS simple ATT close to truth")
check(abs(sa$att - true_att) < 0.15, "SA ATT close to truth")
check(abs(tw$att - true_att) > abs(cs$att - true_att), "TWFE is further from truth than CS (heterogeneity bias)")
for (r in list(tw, cs, sa)) {
    check(all(c(-4:4) %in% r$dynamic$e), paste(r$method, "dynamic covers window"))
    check(r$dynamic$est[r$dynamic$e == -1] == 0, paste(r$method, "reference period is 0"))
}
# true dynamic effect at e=2 averaged over cohorts
true_e2 <- mean(sim$tau[sim$g > 0 & sim$year - sim$g == 2])
est_e2 <- c(cs = cs$dynamic$est[cs$dynamic$e == 2], sa = sa$dynamic$est[sa$dynamic$e == 2])
cat("true e=2:", round(true_e2, 3), " estimates:", round(est_e2, 3), "\n")
check(all(abs(est_e2 - true_e2) < 0.3), "CS and SA dynamic effect at e=2 near truth")

## ---- own SA implementation must match fixest::sunab (plain R only) ----
dfs <- d$data; dfs$dc_gs <- ifelse(dfs$dc_g > 0, dfs$dc_g, 10000)
ref <- fixest::feols(dc_y ~ sunab(dc_gs, dc_t, ref.p = -1) | dc_id + dc_t, data = dfs,
                     vcov = ~dc_id, notes = FALSE, nthreads = 1)
ca <- fixest::coeftable(summary(ref, agg = "att")); ce <- fixest::coeftable(summary(ref))
cat(sprintf("SA own ATT %.6f (se %.6f) vs sunab %.6f (se %.6f)\n", sa$att, sa$se, ca[1, 1], ca[1, 2]))
check(abs(sa$att - ca[1, 1]) < 1e-6, "SA overall ATT equals fixest::sunab")
check(abs(sa$se - ca[1, 2]) < 1e-4, "SA overall SE equals fixest::sunab")
e_ref <- as.numeric(sub(".*::", "", rownames(ce)))
for (k in c(-3, 0, 2, 4)) {
    own <- sa$dynamic$est[sa$dynamic$e == k]; fx <- ce[e_ref == k, 1]
    check(abs(own - fx) < 1e-6, sprintf("SA event time %d equals fixest::sunab", k))
}

## ---- error paths ----
bad <- sim; bad$D[bad$id == 1 & bad$year == 12] <- 0
bad$D[bad$id == which(cohort == 5)[1] & bad$year == 12] <- 0
check(inherits(try(dc_prepare(bad, "y", "id", "year", "D", "indicator"), silent = TRUE), "try-error"),
      "non-absorbing treatment is rejected")
nn <- sim[sim$g > 0, ]
r <- dc_cs(dc_prepare(nn, "y", "id", "year", "D", "indicator"), -4, 4, "nevertreated")
check(!r$ok && nzchar(r$message), "CS without never-treated reports an error instead of crashing")

## ---- M2 estimators: plain R first, then under jamovi's formula sandbox ----
plain <- list(tw = tw, cs = cs, sa = sa,
              bjs = dc_bjs(d, -4, 4), gd = dc_gardner(d, -4, 4), st = dc_stacked(d, -4, 4))
for (r in plain[c("bjs", "gd", "st")])
    cat(sprintf("%-26s ATT %.3f (se %.3f) ok=%s %s\n", r$method, r$att, r$se, r$ok, r$message))
check(all(vapply(plain, function(r) r$ok, TRUE)), "BJS, Gardner, stacked ran")
check(abs(plain$bjs$att - plain$gd$att) < 1e-6, "BJS and Gardner overall ATT agree (both imputation)")
for (k in 0:4) {
    a <- plain$bjs$dynamic$est[plain$bjs$dynamic$e == k]; b <- plain$gd$dynamic$est[plain$gd$dynamic$e == k]
    check(abs(a - b) < 1e-6, sprintf("BJS and Gardner agree at e=%d", k))
}
check(abs(plain$bjs$att - true_att) < 0.2, "BJS ATT near truth")
check(abs(plain$st$att - true_att) < 0.3, "Stacked ATT near truth")
check(all(c(-4:4) %in% plain$st$dynamic$e), "Stacked dynamic covers window")

source(file.path(root, "dev/sandbox.R")); dc_enable_sandbox()
check(inherits(try(stats::as.formula("y ~ sunab(a, b)"), silent = TRUE), "try-error"),
      "sandbox emulation is active")
under <- list(tw = dc_twfe(d, -4, 4), cs = dc_cs(d, -4, 4), sa = dc_sa(d, -4, 4),
              bjs = dc_bjs(d, -4, 4), gd = dc_gardner(d, -4, 4), st = dc_stacked(d, -4, 4))
for (k in names(under)) {
    check(under[[k]]$ok, paste(under[[k]]$method, "runs under sandbox", under[[k]]$message))
    check(abs(under[[k]]$att - plain[[k]]$att) < 1e-8, paste(under[[k]]$method, "sandbox result equals plain R"))
}

## ---- mpdta ----
suppressMessages(library(did)); data(mpdta, package = "did")
dm <- dc_prepare(mpdta, "lemp", "countyreal", "year", "first.treat", "cohort")
csm <- dc_cs(dm, -3, 2)
ref <- aggte(att_gt("lemp", "year", "countyreal", "first.treat", data = mpdta,
                    bstrap = FALSE, cband = FALSE), type = "simple", bstrap = FALSE, cband = FALSE)
cat(sprintf("mpdta CS ATT %.5f vs direct did %.5f\n", csm$att, ref$overall.att))
check(abs(csm$att - ref$overall.att) < 1e-8, "mpdta: adapter reproduces did::aggte simple ATT")
cat("\nAll core tests passed.\n")

## ---- M3 diagnostics (sandbox still active) ----
cat("\n--- diagnostics ---\n")
bc <- dc_bacon(d)
check(bc$ok, paste("Bacon runs", bc$message))
check(abs(sum(bc$detail$weight) - 1) < 1e-8, "Bacon weights sum to 1")
check(abs(bc$twfe - tw$att) < 1e-6, "Bacon-implied TWFE equals the TWFE coefficient")
check(all(c("Treated vs Untreated", "Earlier vs Later Treated", "Later vs Earlier Treated") %in% bc$summary$type),
      "Bacon finds all three comparison types in staggered data")
unbal <- d; unbal$data <- unbal$data[-5, ]; unbal$info$balanced <- FALSE
check(!dc_bacon(unbal)$ok, "Bacon refuses unbalanced panels with a message")

wt <- dc_twfe_weights(d)
check(wt$ok && abs(sum(wt$w) - 1) < 1e-8, "TWFE weights sum to 1 over treated cells")
check(wt$neg_n > 0 && wt$neg_w < 0, "staggered design gives some negative TWFE weights")
check(abs(wt$pos_w + wt$neg_w - 1) < 1e-8, "positive and negative weights add to 1")
# exact identity: with noise-free Y = FE + tau_it * D, beta_TWFE = sum_w tau
dd <- d; set.seed(2)
tau_it <- rnorm(nrow(dd$data), 1, .5)
dd$data$dc_y <- rnorm(N)[dd$data$dc_id] + rnorm(T)[dd$data$dc_t] + tau_it * dd$data$dc_d
beta <- fixest::coeftable(fixest::feols(dc_y ~ dc_d | dc_id + dc_t, data = dd$data, notes = FALSE))[1, 1]
w2 <- dc_twfe_weights(dd)$w
check(abs(beta - sum(w2 * tau_it[dd$data$dc_d == 1])) < 1e-8, "TWFE coefficient = sum of weights x cell effects (exact)")

pt <- dc_pretests(d, -4)
print(pt[, 1:4])
check(all(is.finite(pt$p)) && all(pt$p >= 0 & pt$p <= 1), "pre-tests return p-values for TWFE, SA and CS")
# CS df and p must be internally consistent with the Wald statistic
cs_row <- pt[pt$method == "Callaway-Sant'Anna", ]
check(abs(stats::pchisq(cs_row$stat, cs_row$df, lower.tail = FALSE) - cs_row$p) < 0.02,
      "CS pre-test df is consistent with its reported p-value")
check(pt$p[pt$method == "TWFE"] < 0.01, "TWFE pre-test rejects (TWFE is contaminated by heterogeneity)")
check(pt$p[pt$method == "Sun-Abraham"] > 0.05, "SA pre-test does not reject (parallel trends hold by construction)")
nopre <- dc_pretests(d, -1)
check(all(is.na(nopre$p[nopre$method != "Callaway-Sant'Anna"])) && all(nzchar(nopre$note[nopre$method != "Callaway-Sant'Anna"])),
      "TWFE and SA pre-tests explain themselves when the window has no pre-periods")
check(is.finite(nopre$p[nopre$method == "Callaway-Sant'Anna"]), "CS pre-test uses all pre-periods regardless of the window")

prep <- dc_honest_prepare(d, "sa", -4, 4)
check(prep$npre == 3 && prep$npost == 5, "HonestDiD prep: 3 pre and 5 post periods")
orig <- dc_honest_one(prep, NA, "first", 0.05)
m0 <- dc_honest_one(prep, 0, "first", 0.05); m2 <- dc_honest_one(prep, 2, "first", 0.05)
cat(sprintf("HonestDiD first post: original [%.3f, %.3f]  Mbar=0 [%.3f, %.3f]  Mbar=2 [%.3f, %.3f]\n",
            orig[1], orig[2], m0[1], m0[2], m2[1], m2[2]))
check(m2[["ub"]] - m2[["lb"]] > m0[["ub"]] - m0[["lb"]], "HonestDiD interval widens with Mbar")
check(m0[["lb"]] <= sa$dynamic$est[sa$dynamic$e == 0] && m0[["ub"]] >= sa$dynamic$est[sa$dynamic$e == 0],
      "HonestDiD interval contains the point estimate")
prep_cs <- dc_honest_prepare(d, "cs", -4, 4)
check(prep_cs$npre == 3 && prep_cs$npost == 5, "HonestDiD prep via CS influence functions works")
check(inherits(try(dc_honest_prepare(d, "sa", -1, 4), silent = TRUE), "try-error"),
      "HonestDiD refuses a window without extra pre-periods")
cat("\nAll diagnostics tests passed.\n")

## ---- M4 covariates (sandbox still active) ----
cat("\n--- covariates ---\n")
set.seed(5); N2 <- 1000; T2 <- 8
x <- rnorm(N2); grp <- sample(c("a", "b", "c"), N2, TRUE)
g2 <- ifelse(runif(N2) < .65, ifelse(runif(N2) < plogis(1.5 * x), 4, 6), 0)   # cohort depends on x
cv <- expand.grid(id = 1:N2, year = 1:T2)
cv$g <- g2[cv$id]; cv$x <- x[cv$id]; cv$grp <- grp[cv$id]
cv$z <- rnorm(nrow(cv))                                  # time-varying covariate
cv$D <- as.integer(cv$g > 0 & cv$year >= cv$g)
cv$y <- rnorm(N2)[cv$id] + rnorm(T2, sd = .3)[cv$year] + 0.4 * cv$x * cv$year + 0.5 * cv$z +
        1 * cv$D + rnorm(nrow(cv), sd = .5)               # true effect = 1
d0 <- dc_prepare(cv, "y", "id", "year", "D", "indicator")
dx <- dc_prepare(cv, "y", "id", "year", "D", "indicator", covs = c("x", "z"))
check(nrow(dx$info$covs) == 2 && all(c("dc_x1", "dc_x2") %in% names(dx$data)), "prepare: numeric covariates mapped to dc_x1, dc_x2")
dg <- dc_prepare(cv, "y", "id", "year", "D", "indicator", covs = "grp")
check(nrow(dg$info$covs) == 2 && grepl("grp = b", dg$info$covs$label[1]), "prepare: factor covariate becomes dummies")
check(grepl("Covariates", dc_info_html(dx$info)), "info HTML lists covariates")
cs_no <- dc_cs(d0, -3, 3); cs_x <- dc_cs(dx, -3, 3, "nevertreated", "dr")
cat(sprintf("covariate DGP (true 1): CS no covs %.3f, CS with x,z (dr) %.3f\n", cs_no$att, cs_x$att))
check(cs_no$ok && cs_x$ok, "CS runs with and without covariates")
check(abs(cs_x$att - 1) < abs(cs_no$att - 1), "CS covariate adjustment moves the estimate toward the truth")
check(abs(cs_x$att - 1) < 0.2, "CS with covariates is close to the true effect")
for (m in c("ipw", "reg")) check(dc_cs(dx, -3, 3, "nevertreated", m)$ok, paste("CS estimation method", m, "runs"))
allx <- list(dc_twfe(dx, -3, 3), dc_sa(dx, -3, 3), dc_bjs(dx, -3, 3), dc_gardner(dx, -3, 3), dc_stacked(dx, -3, 3))
for (r in allx) check(isTRUE(r$ok) && is.finite(r$att) && all((if (r$method == "Borusyak-Jaravel-Spiess") c(-1, 0:3) else -3:3) %in% r$dynamic$e),
                      paste(r$method, "runs with covariates under the sandbox", r$message))
ally <- list(dc_twfe(dg, -3, 3), dc_sa(dg, -3, 3), dc_bjs(dg, -3, 3), dc_gardner(dg, -3, 3), dc_stacked(dg, -3, 3), dc_cs(dg, -3, 3))
for (r in ally) check(r$ok, paste(r$method, "runs with dummy-coded factor covariate", r$message))
# the time-varying covariate z (true coefficient .5) should change TWFE-type estimates only slightly but must be used
tw_nz <- dc_twfe(d0, -3, 3)$att; tw_z <- allx[[1]]$att
cat(sprintf("TWFE without covariates %.3f, with x,z %.3f\n", tw_nz, tw_z))
miss <- cv; miss$z[sample(nrow(miss), 50)] <- NA
dm <- dc_prepare(miss, "y", "id", "year", "D", "indicator", covs = c("x", "z"))
check(dm$info$n_missing == 50, "rows with missing covariate values are dropped and counted")
cat("\nAll covariate tests passed.\n")

## ---- M4 triple differences (sandbox still active) ----
cat("\n--- DDD ---\n")
dd3 <- read.csv(file.path(root, "dev/csv/sim_ddd.csv"))
true_ddd <- 1.489488   # mean effect on treated eligible observations (see data generation)
dp <- dc_prepare(dd3, "y", "id", "year", "g", "cohort", part = "elig")
check("dc_p" %in% names(dp$data) && all(dp$data$dc_p %in% 0:1), "prepare: partition variable carried through")
ddd <- dc_ddd(dp, -3, 3, "nevertreated", "dr")
el <- dc_did_part(dp, "eligible", -3, 3, "nevertreated", "dr")
inel <- dc_did_part(dp, "ineligible", -3, 3, "nevertreated", "dr")
for (r in list(ddd, el, inel))
    cat(sprintf("%-26s ATT %.3f (se %.3f) ok=%s %s\n", r$method, r$att, r$se, r$ok, r$message))
check(ddd$ok && el$ok && inel$ok, "DDD and the two DiD comparisons run under the sandbox")
check(abs(ddd$att - true_ddd) < 0.2, "DDD recovers the effect on eligible units")
check(el$att - true_ddd > 0.5, "DiD on eligible units alone is biased upward by the group shock")
check(abs(inel$att - 0.8) < 0.15, "DiD on ineligible units picks up the 0.8 group-wide shock")
check(all(-3:3 %in% ddd$dynamic$e), "DDD event study covers the window")
check(abs(ddd$dynamic$est[ddd$dynamic$e == -1]) < 1e-12, "DDD reference period is 0")
pre <- ddd$dynamic[ddd$dynamic$e < -1, ]
check(all(abs(pre$est) < 3 * pre$se), "DDD pre-period effects are close to zero")
bad <- dd3; bad$elig[1] <- 2
check(inherits(try(dc_prepare(bad, "y", "id", "year", "g", "cohort", part = "elig"), silent = TRUE), "try-error"),
      "partition values other than 0/1 are rejected")
bad2 <- dd3; bad2$elig[bad2$id == 1 & bad2$year == 3] <- 1 - bad2$elig[bad2$id == 1 & bad2$year == 3]
check(inherits(try(dc_prepare(bad2, "y", "id", "year", "g", "cohort", part = "elig"), silent = TRUE), "try-error"),
      "partition that varies within a unit is rejected")
cat("\nAll DDD tests passed.\n")

## ---- M5 cluster, bootstrap, unbalanced, BJS pre-periods (sandbox still active) ----
cat("\n--- M5 core ---\n")
# BJS pre-period coefficients on the main simulated panel (parallel trends hold)
bjs5 <- dc_bjs(d, -4, 4)
check(bjs5$ok && all(-4:4 %in% bjs5$dynamic$e), "BJS now reports pre-treatment coefficients")
pre5 <- bjs5$dynamic[bjs5$dynamic$e < -1, ]
check(all(abs(pre5$est) < 3 * pre5$se), "BJS pre-period coefficients are close to zero under parallel trends")
check(abs(bjs5$att - plain$bjs$att) < 1e-10, "BJS overall ATT unchanged by the pre-period addition")

# cluster variable
set.seed(11)
sim$cl <- (sim$id - 1) %/% 10 + 1        # clusters of 10 units
dcl <- dc_prepare(sim, "y", "id", "year", "D", "indicator", cluster = "cl")
check(isTRUE(dcl$info$cluster) && length(unique(dcl$data$dc_cl)) == N / 10, "prepare: cluster variable mapped")
badcl <- sim; badcl$cl[badcl$id == 1 & badcl$year == 2] <- 999
check(inherits(try(dc_prepare(badcl, "y", "id", "year", "D", "indicator", cluster = "cl"), silent = TRUE), "try-error"),
      "cluster variable that varies within a unit is rejected")
tw_u <- dc_twfe(d, -4, 4); tw_c <- dc_twfe(dcl, -4, 4)
check(abs(tw_u$att - tw_c$att) < 1e-10 && tw_c$se != tw_u$se, "clustering changes TWFE standard errors but not the estimate")
for (fn in list(function(x) dc_sa(x, -4, 4), function(x) dc_bjs(x, -4, 4),
                function(x) dc_gardner(x, -4, 4), function(x) dc_stacked(x, -4, 4))) {
    a <- fn(d); b <- fn(dcl)
    check(a$ok && b$ok && abs(a$att - b$att) < 1e-8, paste(a$method, "runs with a cluster variable, same estimate"))
}

# bootstrap for CS: reproducible with a seed, close to analytic SE, forced by clustering
cs_a <- dc_cs(d, -4, 4)
b1 <- dc_cs(d, -4, 4, boot = list(on = TRUE, reps = 300, seed = 7))
b2 <- dc_cs(d, -4, 4, boot = list(on = TRUE, reps = 300, seed = 7))
check(b1$ok && b1$att == cs_a$att, "CS bootstrap leaves the point estimate unchanged")
check(identical(b1$se, b2$se), "CS bootstrap is reproducible for a given seed")
check(abs(b1$se / cs_a$se - 1) < 0.35, sprintf("CS bootstrap SE (%.3f) is close to the analytic SE (%.3f)", b1$se, cs_a$se))
ccl <- dc_cs(dcl, -4, 4, boot = list(on = FALSE, reps = 200, seed = 7))
check(ccl$ok && grepl("bootstrap", ccl$control), "clustering above the unit level switches CS to the bootstrap")

# unbalanced panel: did's unbalanced-panel path can crash R, so CS uses the complete units
set.seed(12)
att_units <- sample(unique(sim$id), 50)                    # 50 units lose their last 1-3 periods
ub <- sim[!(sim$id %in% att_units & sim$year > T - sample(1:3, nrow(sim), TRUE)), ]
dub <- dc_prepare(ub, "y", "id", "year", "D", "indicator")
check(!dub$info$balanced, "unbalanced panel is detected")
csub <- dc_cs(dub, -4, 4)
check(csub$ok && grepl("set aside", csub$control), "CS runs on the complete units and says how many were set aside")
check(abs(csub$att - true_att) < 0.25, sprintf("CS on a mildly unbalanced panel is close to the truth (ATT %.3f)", csub$att))
heavy <- sim[-sample(nrow(sim), 1500), ]                   # almost every unit misses some period
dh <- dc_prepare(heavy, "y", "id", "year", "D", "indicator")
rh <- dc_cs(dh, -4, 4)
check(inherits(rh, "list") && (rh$ok || nzchar(rh$message)), "heavily unbalanced panel returns a result or a message, never a crash")
check(all(c(dc_twfe(dh, -4, 4)$ok, dc_sa(dh, -4, 4)$ok, dc_bjs(dh, -4, 4)$ok)), "fixed-effects estimators still use the unbalanced panel")

# CS extras
cse <- dc_cs(d, -4, 4, extras = TRUE)
check(!is.null(cse$extra$gt) && nrow(cse$extra$gt) > 10 && all(c("g", "t", "att", "se") %in% names(cse$extra$gt)),
      "CS returns group-time effects")
check(setequal(cse$extra$group$g, c(5, 7, 9)), "CS group aggregation covers all cohorts")
check(nrow(cse$extra$calendar) >= 5, "CS calendar aggregation returned")
cat("\nAll M5 core tests passed.\n")

## ---- M6 time coding and input validation (sandbox still active) ----
cat("\n--- M6 time coding and validation ---\n")
base_res <- list(plain$tw, plain$cs, plain$sa, plain$bjs, plain$gd, plain$st)
att_of <- function(dd) vapply(list(dc_twfe(dd, -4, 4), dc_cs(dd, -4, 4), dc_sa(dd, -4, 4), dc_bjs(dd, -4, 4),
                                   dc_gardner(dd, -4, 4), dc_stacked(dd, -4, 4)), function(r) if (isTRUE(r$ok)) r$att else NA_real_, 0)
base_att <- vapply(base_res, function(r) r$att, 0)
variants <- list(
    from_zero  = transform(sim, year = year - 1),
    negative   = transform(sim, year = year - 6),
    step2      = transform(sim, year = 2000 + 2 * (year - 1)),
    step5_neg  = transform(sim, year = -20 + 5 * (year - 1)))
for (nm in names(variants)) {
    dv <- dc_prepare(variants[[nm]], "y", "id", "year", "D", "indicator")
    a <- att_of(dv)
    check(all(is.finite(a)) && max(abs(a - base_att)) < 1e-8, sprintf("time coding '%s': all six estimators equal the 1..T coding", nm))
    check(isTRUE(dv$info$time_map$shifted), sprintf("time coding '%s' is recoded", nm))
}
v2 <- dc_prepare(variants$step2, "y", "id", "year", "D", "indicator")
html <- dc_info_html(v2$info)
check(grepl("2000 to 2022", html) && grepl("apart", html), "info shows the original periods and explains the recoding (step 2)")
check(setequal(dc_orig_time(v2$info, as.numeric(names(v2$info$cohorts))), c(2008, 2012, 2016)), "cohorts are shown in original periods (2008, 2012, 2016)")
vn <- dc_prepare(variants$negative, "y", "id", "year", "D", "indicator")
check(grepl("-5 to 6", dc_info_html(vn$info)) && grepl("start at -5", dc_info_html(vn$info)), "negative periods: original values shown with an explanation")

# first treated period given as a variable, original period units
sc <- transform(variants$step2, g2 = ifelse(g > 0, 2000 + 2 * (g - 1), 0))
dsc <- dc_prepare(sc, "y", "id", "year", "g2", "cohort")
check(max(abs(att_of(dsc) - base_att)) < 1e-8, "cohort input in original periods (step 2) equals the standard coding")
off <- transform(sc, g2 = ifelse(g2 > 0, g2 + 1, 0))         # odd years do not exist in the panel
check(inherits(try(dc_prepare(off, "y", "id", "year", "g2", "cohort"), silent = TRUE), "try-error"), "first treated period that is not an observed period is rejected")
irr <- sim; irr$year <- c(2000, 2001, 2003, 2007, 2008, 2012, 2013, 2020, 2021, 2022, 2025, 2026)[irr$year]
check(isTRUE(dc_prepare(irr, "y", "id", "year", "D", "indicator")$info$time_map$irregular), "irregular spacing is flagged")

# a unit first treated at period 0 used to be misread as never treated
z <- transform(sim, year = year - 5)                          # periods -4..7: cohort 5 starts at 0
dz <- dc_prepare(z, "y", "id", "year", "D", "indicator")
check(dz$info$n_never == d$info$n_never && sum(dz$info$cohorts) == sum(d$info$cohorts) && dz$info$n_always == 0,
      "units first treated at period 0 are not mistaken for never-treated")

# friendly messages
msg <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
tt <- sim; tt$D <- ifelse(tt$D == 1, "yes", "no")
check(grepl("numeric 0/1", msg(dc_prepare(tt, "y", "id", "year", "D", "indicator"))), "text treatment indicator: clear message")
td <- sim; td$year <- sprintf("%d-01-01", 1999 + td$year)
check(grepl("Dates and text", msg(dc_prepare(td, "y", "id", "year", "D", "indicator"))), "date time variable: clear message")
nt <- sim; nt$D <- 0L
check(grepl("No unit is treated", msg(dc_prepare(nt, "y", "id", "year", "D", "indicator"))), "no treated units: clear message")
cy <- sim; cy$y <- 3
check(grepl("constant", msg(dc_prepare(cy, "y", "id", "year", "D", "indicator"))), "constant outcome: clear message")
cat("\nAll M6 core tests passed.\n")

## ---- M6 fast Goodman-Bacon decomposition equals bacondecomp::bacon exactly ----
cat("\n--- fast Bacon ---\n")
cmp_bacon <- function(dd, label) {
    slow <- as.data.frame(bacondecomp::bacon(dc_y ~ dc_d, data = dd$data, id_var = "dc_id", time_var = "dc_t", quietly = TRUE))
    fast <- dc_bacon_2x2(dd$data)
    key <- function(x) paste(x$treated, x$untreated)
    slow <- slow[order(key(slow)), ]; fast <- fast[order(key(fast)), ]
    check(identical(key(slow), key(fast)) && all(slow$type == fast$type), paste(label, ": same comparisons and types"))
    check(max(abs(slow$weight - fast$weight)) < 1e-10, paste(label, ": weights identical"))
    check(max(abs(slow$estimate - fast$estimate)) < 1e-8, paste(label, ": 2x2 estimates identical"))
}
cmp_bacon(d, "simulated staggered panel")
cmp_bacon(dc_prepare(as.data.frame(bacondecomp::castle), "l_homicide", "state", "year", "post", "indicator"), "Castle Doctrine")
set.seed(21)
rg <- expand.grid(id = 1:90, t = 1:9); gg <- sample(c(2, 4, 5, 8, 0), 90, TRUE, prob = c(.1, .3, .2, .15, .25))
rg$g <- gg[rg$id]; rg$D <- as.integer(rg$g > 0 & rg$t >= rg$g); rg$y <- rnorm(90)[rg$id] + .3 * rg$t + rg$D * (1 + .2 * (rg$t - rg$g)) + rnorm(nrow(rg))
cmp_bacon(dc_prepare(rg, "y", "id", "t", "D", "indicator"), "random design with unequal cohort sizes")
nn <- rg[rg$g != 0, ]
cmp_bacon(dc_prepare(nn, "y", "id", "t", "D", "indicator"), "design without never-treated units")
cat("\nAll fast-Bacon tests passed.\n")

## ---- M7 economists' standards: t(G-1) inference, sample sizes, sampling weights (sandbox still active) ----
cat("\n--- M7 inference, N, weights ---\n")
# (1) regression-based estimators report t(G-1) inference exactly as fixest does
set.seed(3); Ns <- 12; Ts <- 10
gs <- sample(c(5, 7, 0), Ns, TRUE); smx <- expand.grid(id = 1:Ns, year = 1:Ts); smx$g <- gs[smx$id]
smx$D <- as.integer(smx$g > 0 & smx$year >= smx$g); smx$y <- rnorm(Ns)[smx$id] + .1 * smx$year + .5 * smx$D + rnorm(nrow(smx), sd = .6)
ds <- dc_prepare(smx, "y", "id", "year", "D", "indicator")
tws <- dc_twfe(ds, -3, 3)
fx <- fixest::coeftable(fixest::feols(dc_y ~ dc_d | dc_id + dc_t, data = ds$data, vcov = ~dc_cl, notes = FALSE))
inf <- dc_inference(tws$att, tws$se, tws$df, 95)
check(tws$df == Ns - 1, "TWFE reports t(G-1) degrees of freedom (G = 12 clusters)")
check(abs(inf$p - fx[1, 4]) < 1e-12, sprintf("TWFE p-value equals fixest's (%.5f); the normal approximation would give %.5f", inf$p, 2 * pnorm(-abs(tws$att / tws$se))))
check(inf$crit > stats::qnorm(.975), "t critical value is wider than the normal one with 12 clusters")
check(is.na(dc_cs(ds, -3, 3)$df) && is.na(dc_bjs(ds, -3, 3)$df), "CS and BJS keep their packages' normal inference")
check(abs(dc_inference(1, .1, NA, 95)$crit - qnorm(.975)) < 1e-12, "dc_inference falls back to the normal distribution when df is NA")

# (2) sample sizes
r6 <- list(dc_twfe(d, -4, 4), dc_cs(d, -4, 4), dc_sa(d, -4, 4), dc_bjs(d, -4, 4), dc_gardner(d, -4, 4), dc_stacked(d, -4, 4))
check(all(vapply(r6, function(r) is.finite(r$n_obs) && r$n_obs > 0, TRUE)), "every estimator reports its number of observations")
check(r6[[1]]$n_obs == nrow(d$data) && r6[[1]]$n_cl == N, "TWFE: N observations and G clusters")
check(r6[[6]]$n_obs > r6[[1]]$n_obs, "stacked DiD reports the larger stacked sample")

# (3) sampling weights: a weight of w must equal w copies of the unit (frequency-weight equivalence)
set.seed(8); Nw <- 60; Tw <- 8
gw <- sample(c(4, 6, 0), Nw, TRUE, prob = c(.3, .3, .4)); ww <- sample(1:3, Nw, TRUE)
wp <- expand.grid(id = 1:Nw, year = 1:Tw); wp$g <- gw[wp$id]; wp$w <- ww[wp$id]
wp$D <- as.integer(wp$g > 0 & wp$year >= wp$g)
wp$y <- rnorm(Nw)[wp$id] + .2 * wp$year + wp$D * (1 + .2 * (wp$year - wp$g)) + rnorm(nrow(wp), sd = .5)
rep_ids <- unlist(lapply(1:Nw, function(i) rep(i, ww[i])))
rp <- do.call(rbind, lapply(seq_along(rep_ids), function(k) { x <- wp[wp$id == rep_ids[k], ]; x$id <- k; x }))
dW <- dc_prepare(wp, "y", "id", "year", "D", "indicator", weights = "w")
dR <- dc_prepare(rp, "y", "id", "year", "D", "indicator")
check(isTRUE(dW$info$weighted) && all(dR$data$dc_sw == 1), "prepare: weights flagged, unweighted data carry weight 1")
for (nm in c("twfe", "sa", "cs", "gardner", "stacked")) {
    fn <- switch(nm, twfe = dc_twfe, sa = dc_sa, cs = dc_cs, gardner = dc_gardner, stacked = dc_stacked)
    a <- fn(dW, -3, 3); b <- fn(dR, -3, 3)
    check(isTRUE(a$ok) && isTRUE(b$ok), paste(nm, "runs with and without weights", a$message, b$message))
    check(abs(a$att - b$att) < 1e-7, sprintf("%s: weighted ATT equals the replicated-unit ATT (%.6f)", nm, a$att))
    ea <- a$dynamic; eb <- b$dynamic; k <- match(ea$e, eb$e)
    check(max(abs(ea$est - eb$est[k])) < 1e-7, paste(nm, ": weighted event-time estimates equal the replicated-unit ones"))
}
bw <- dc_bjs(dW, -3, 3)
check(!bw$ok && grepl("sampling weights", bw$message), "BJS explains that it cannot use sampling weights")
unw <- dc_twfe(dc_prepare(wp, "y", "id", "year", "D", "indicator"), -3, 3)$att
check(abs(unw - dc_twfe(dW, -3, 3)$att) > 1e-4, "weights actually change the TWFE estimate")
w1 <- wp; w1$w <- 2.5
d1 <- dc_prepare(wp, "y", "id", "year", "D", "indicator"); d1w <- dc_prepare(w1, "y", "id", "year", "D", "indicator", weights = "w")
check(abs(dc_twfe(d1, -3, 3)$att - dc_twfe(d1w, -3, 3)$att) < 1e-9 && abs(dc_sa(d1, -3, 3)$att - dc_sa(d1w, -3, 3)$att) < 1e-9, "constant weights change nothing")
msg2 <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
nw <- wp; nw$w <- -1
check(grepl("positive", msg2(dc_prepare(nw, "y", "id", "year", "D", "indicator", weights = "w"))), "non-positive weights are rejected")
tv <- wp; tv$w <- tv$year
check(grepl("constant within", msg2(dc_prepare(tv, "y", "id", "year", "D", "indicator", weights = "w"))), "time-varying weights are rejected")
tw2 <- wp; tw2$w <- "heavy"
check(grepl("numeric", msg2(dc_prepare(tw2, "y", "id", "year", "D", "indicator", weights = "w"))), "text weights are rejected")

# (4) small-sample warnings in the panel summary
html_small <- dc_info_html(ds$info)
check(grepl("Only 12 clusters", html_small), "few clusters: the summary warns")
check(!grepl("Only", dc_info_html(d$info)) , "no few-cluster warning with 600 clusters")
check(grepl("Sampling weights: w", dc_info_html(dW$info)), "the summary names the weight variable")
cat("\nAll M7 tests passed.\n")

## ---- M7 the bootstrap must not disturb the caller's random-number state ----
set.seed(99); ref_draw <- runif(3)
set.seed(99)
bb <- dc_cs(d, -4, 4, boot = list(on = TRUE, reps = 99, seed = 5))
after <- runif(3)
check(isTRUE(bb$ok) && identical(ref_draw, after), "CS bootstrap leaves the caller's RNG state untouched")
rm(.Random.seed, envir = globalenv())
bb2 <- dc_cs(d, -4, 4, boot = list(on = TRUE, reps = 99, seed = 5))
check(isTRUE(bb2$ok) && !exists(".Random.seed", envir = globalenv(), inherits = FALSE), "no RNG state is created when the caller had none")
cat("\nAll RNG-state tests passed.\n")

## ---- M7 infinite outcomes are dropped like missing values ----
inf_d <- sim; inf_d$y[c(5, 77)] <- c(Inf, -Inf)
di <- dc_prepare(inf_d, "y", "id", "year", "D", "indicator")
check(di$info$n_missing == 2 && all(is.finite(di$data$dc_y)), "infinite outcomes are dropped and counted")
check(isTRUE(dc_twfe(di, -4, 4)$ok), "estimation proceeds after dropping infinite outcomes")
cat("\nAll infinite-outcome tests passed.\n")
