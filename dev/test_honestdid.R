# Vendored HonestDiD subset vs the HonestDiD package (identical bounds).
root <- if (file.exists("R/core_prepare.R")) "." else ".."
source(file.path(root, "R/vendor_honestdid.R"))
check <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) quit(status = 1) }
has_ref <- requireNamespace("HonestDiD", quietly = TRUE)
set.seed(11)
for (case in list(c(4, 3), c(3, 2), c(2, 4))) {
    npre <- case[1]; npost <- case[2]; k <- npre + npost
    A <- matrix(rnorm(k * k), k); sigma <- crossprod(A) / 20 + diag(k) * .01
    beta <- c(rnorm(npre, 0, .05), 0.4 + cumsum(rnorm(npost, .1, .05)))
    for (l in list(hd_basisVector(1, npost), rep(1 / npost, npost))) {
        a <- hd_constructOriginalCS(beta, sigma, npre, npost, l_vec = l, alpha = .05)
        b <- hd_createSensitivityResults_relativeMagnitudes(beta, sigma, numPrePeriods = npre, numPostPeriods = npost,
                                                            Mbarvec = c(0.5, 1, 2), l_vec = l, alpha = .05)
        check(all(is.finite(c(a$lb, a$ub, b$lb, b$ub))), sprintf("finite bounds (npre %d, npost %d)", npre, npost))
        if (has_ref) {
            a2 <- HonestDiD::constructOriginalCS(beta, sigma, npre, npost, l_vec = l, alpha = .05)
            b2 <- HonestDiD::createSensitivityResults_relativeMagnitudes(beta, sigma, numPrePeriods = npre, numPostPeriods = npost,
                                                                         Mbarvec = c(0.5, 1, 2), l_vec = l, alpha = .05)
            check(isTRUE(all.equal(c(a$lb, a$ub), c(a2$lb, a2$ub), tolerance = 1e-8)), "original CS equals HonestDiD")
            check(isTRUE(all.equal(c(b$lb, b$ub), c(b2$lb, b2$ub), tolerance = 1e-6)), "relative-magnitude bounds equal HonestDiD")
        }
    }
}
cat("ALL OK\n")
