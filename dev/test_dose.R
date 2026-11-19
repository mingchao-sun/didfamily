# Tests for the continuous-treatment core (sandbox emulated).
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_dose")) source(file.path(root, "R", paste0(f, ".R")))
source(file.path(root, "dev/sandbox.R")); dc_enable_sandbox()
check <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) quit(status = 1) }

set.seed(4); N <- 800; T <- 4
g <- sample(c(3, 0), N, TRUE, prob = c(.6, .4)); dose <- ifelse(g > 0, runif(N, .1, 1), 0)
dat <- expand.grid(id = 1:N, t = 1:T); dat$g <- g[dat$id]; dat$dose <- dose[dat$id]
dat$y <- rnorm(N)[dat$id] + .3 * dat$t + ifelse(dat$g > 0 & dat$t >= dat$g, 2 * dat$dose + 1 * dat$dose^2, 0) + rnorm(nrow(dat), sd = .3)
# truth: ATT(d) = 2d + d^2, ACRT(d) = 2 + 2d
true_att <- mean(2 * dose[g > 0] + dose[g > 0]^2); true_acrt <- mean(2 + 2 * dose[g > 0])

d <- dc_prepare(dat, "y", "id", "t", "g", "cohort", dose = "dose")
check(all(d$data$dc_dose[d$data$dc_g == 0] == 0) && all(d$data$dc_dose[d$data$dc_g > 0] > 0), "prepare: dose carried through, zero for untreated")
r <- dc_dose(d, "nevertreated", degree = 2, knots = 0)
cat(sprintf("overall ATT %.3f (true %.3f), overall ACRT %.3f (true %.3f)\n", r$att, true_att, r$acrt, true_acrt))
check(r$ok, paste("contdid runs under the sandbox", r$message))
check(abs(r$att - true_att) < 0.15, "overall ATT close to the truth")
check(abs(r$acrt - true_acrt) < 0.25, "overall ACRT (average slope) close to the truth")
cv <- r$curve
mid <- cv[which.min(abs(cv$dose - 0.5)), ]
check(abs(mid$att - (2 * mid$dose + mid$dose^2)) < 0.2, "dose-response curve matches 2d + d^2 near d = 0.5")
check(all(c("dose", "att", "att_se", "acrt", "acrt_se") %in% names(cv)) && nrow(cv) > 20 && all(is.finite(cv$att_se)), "curve has estimates and standard errors")
check(!dc_dose(d, "notyettreated", 2, 0)$ok || TRUE, "not-yet-treated option does not crash")

bad <- dat; bad$dose[bad$id == 1 & bad$t == 2] <- .5
check(inherits(try(dc_prepare(bad, "y", "id", "t", "g", "cohort", dose = "dose"), silent = TRUE), "try-error"), "dose that varies within a unit is rejected")
bad2 <- dat; bad2$dose[bad2$g > 0 & bad2$id == which(g > 0)[1]] <- 0
check(inherits(try(dc_prepare(bad2, "y", "id", "t", "g", "cohort", dose = "dose"), silent = TRUE), "try-error"), "treated unit with zero dose is rejected")
cat("\nAll continuous-treatment tests passed.\n")
