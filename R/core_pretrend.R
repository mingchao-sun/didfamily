# Core layer: power of the conventional pre-trend test (Roth, 2022, "Pretest with caution").
# The pre-test passes when no pre-treatment event-study coefficient is individually significant (the "no individually significant" test);
# it rejects if any |beta_k| / se_k exceeds the critical value. Power against a hypothesized violation of parallel trends is the probability of
# rejecting when the pre-period coefficients equal that violation. Implemented from the paper and the definitions of the R package `pretrends`
# (power from a multivariate normal box probability; the violation is a linear trend through the reference period).

# Probability that the pre-test rejects when the pre-period coefficients have mean `delta` and covariance `sigma` (mvtnorm box probability).
dc_pt_reject <- function(delta, sigma, crit = stats::qnorm(0.975)) {
    sigma <- as.matrix(sigma)
    ub <- sqrt(diag(sigma)) * crit
    if (length(ub) == 1) return(1 - (stats::pnorm(ub, delta, sqrt(sigma[1, 1])) - stats::pnorm(-ub, delta, sqrt(sigma[1, 1]))))
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(20261118)                                    # the box probability is computed by a randomised algorithm: fix the seed
    p <- mvtnorm::pmvnorm(lower = -ub, upper = ub, mean = delta, sigma = sigma, algorithm = mvtnorm::GenzBretz(abseps = 1e-8, maxpts = 250000))
    1 - as.numeric(p)
}

# Slope (per event period) of a linear trend through the reference period that the pre-test rejects with probability `target`.
dc_pt_slope <- function(sigma_pre, e_pre, ref = -1, target = 0.5, crit = stats::qnorm(0.975)) {
    rel <- e_pre - ref
    f <- function(s) dc_pt_reject(rel * s, sigma_pre, crit) - target
    hi <- 8 * max(sqrt(diag(as.matrix(sigma_pre))))
    stats::uniroot(f, c(0, hi), tol = 1e-7)$root
}

# prep: output of dc_honest_prepare (beta, sigma, e without the reference period). slope: hypothesized trend slope per period (NA or 0: only the
# slopes for the target powers are computed). Returns the observed pre-test, the power against the slope, the size, the slopes for target powers
# and the bias that each slope implies for the post-treatment periods.
dc_pretrend <- function(prep, slope = NA_real_, level = 0.95, targets = c(0.5, 0.8), ref = -1) {
    crit <- stats::qnorm(1 - (1 - level) / 2)
    pre <- which(prep$e < ref); post <- which(prep$e > ref)
    if (!length(pre)) stop("The pre-trend power needs at least one pre-treatment period besides the reference period.", call. = FALSE)
    sp <- as.matrix(prep$sigma)[pre, pre, drop = FALSE]; bp <- prep$beta[pre]; ep <- prep$e[pre]
    tstat <- bp / sqrt(diag(sp))
    size <- dc_pt_reject(rep(0, length(pre)), sp, crit)
    sl <- vapply(targets, function(t) dc_pt_slope(sp, ep, ref, t, crit), numeric(1))
    power <- if (is.finite(slope) && slope != 0) dc_pt_reject((ep - ref) * slope, sp, crit) else NA_real_
    wald <- as.numeric(t(bp) %*% solve(sp, bp)); wp <- stats::pchisq(wald, length(bp), lower.tail = FALSE)
    list(e_pre = ep, beta_pre = bp, se_pre = sqrt(diag(sp)), t_pre = tstat, passes = all(abs(tstat) <= crit), max_t = max(abs(tstat)), crit = crit,
         size = size, power = power, slope = slope, targets = targets, slope_for_power = sl, wald = wald, wald_p = wp, wald_df = length(bp),
         e_post = prep$e[post], beta_post = prep$beta[post], se_post = sqrt(diag(as.matrix(prep$sigma))[post]), ref = ref)
}

dc_fmt3 <- function(x) formatC(x, format = "g", digits = 3)
