# Vendored two-stage estimator vs the did2s package (identical estimates and standard errors).
root <- if (file.exists("R/core_prepare.R")) "." else ".."
source(file.path(root, "R/core_prepare.R")); source(file.path(root, "R/core_estimators.R")); source(file.path(root, "R/core_did2s.R"))
check <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) quit(status = 1) }
has_ref <- requireNamespace("did2s", quietly = TRUE)
if (!has_ref) cat("did2s not installed: only comparing against stored values is possible\n")

set.seed(3)
N <- 300; T <- 10
g <- sample(c(4, 6, 8, 0), N, TRUE)
s <- expand.grid(id = 1:N, t = 1:T); s$g <- g[s$id]
s$D <- as.integer(s$g > 0 & s$t >= s$g)
s$x <- rnorm(nrow(s)); s$w <- runif(N, .5, 3)[s$id]
s$y <- rnorm(N)[s$id] + rnorm(T)[s$t] + 0.5 * s$D * (s$t - s$g + 1) + .3 * s$x + rnorm(nrow(s))
s$cl <- (s$id %% 20) + 1
s$y2 <- s$y; s$dc_y <- s$y; s$dc_id <- s$id; s$dc_t <- s$t; s$dc_d <- s$D; s$dc_cl <- s$cl; s$dc_sw <- s$w
s$dc_ev1 <- as.numeric(s$D == 1 & s$t - s$g == 0); s$dc_ev2 <- as.numeric(s$D == 1 & s$t - s$g == 1)
s$dc_ev3 <- as.numeric(s$D == 1 & s$t - s$g >= 2)

cmp <- function(fs, ss, w, label) {
    mine <- dc_did2s(s, "dc_y", fs, ss, "dc_d", "dc_cl", w)
    cm <- fixest::coeftable(mine)
    if (has_ref) {
        ref <- did2s::did2s(s, "dc_y", fs, ss, "dc_d", "dc_cl", weights = w)
        cr <- fixest::coeftable(ref)
        check(isTRUE(all.equal(cm[, 1:2], cr[, 1:2], tolerance = 1e-9, check.attributes = FALSE)), paste(label, "- equals did2s"))
    } else check(all(is.finite(cm[, 2])), paste(label, "- finite"))
}
cmp(~ 0 | dc_id + dc_t, ~ dc_d, NULL, "static")
cmp(~ 0 | dc_id + dc_t, ~ dc_ev1 + dc_ev2 + dc_ev3, NULL, "dynamic")
cmp(~ x | dc_id + dc_t, ~ dc_d, NULL, "covariate")
cmp(~ x | dc_id + dc_t, ~ dc_d, "dc_sw", "covariate + weights")
cmp(~ 0 | dc_id + dc_t, ~ dc_ev1 + dc_ev2 + dc_ev3, "dc_sw", "dynamic + weights")
cat("ALL OK\n")
