# Two-stage difference-in-differences (Gardner 2022) with the GMM / influence-function
# variance of Butts & Gardner. This is a port of the estimator in the 'did2s' package
# (Kyle Butts, MIT licence, Copyright 2021 did2s authors): functions did2s(), did2s_estimate()
# and robust_solve_XtX(), analytic-variance branch only. It is vendored so that the module
# does not have to install did2s, whose dependency HonestDiD pulls in CVXR and several
# compiled solvers (gmp, Rglpk, clarabel, highs, osqp, scs). dev/test_did2s.R checks that
# results are identical to did2s.

dc_robust_solve <- function(X, Y) {
    XtY <- if (is.vector(Y)) Matrix::crossprod(X, Y)
           else if (nrow(Y) == nrow(X)) Matrix::crossprod(X, Y)
           else if (nrow(Y) == ncol(X)) Y
           else stop("Incompatible dimensions between X and Y")
    XtX <- Matrix::crossprod(X)
    tryCatch(Matrix::solve(XtX, XtY), error = function(e) {
        SVD <- svd(as.matrix(XtX))
        tol <- max(dim(XtX)) * .Machine$double.eps * max(SVD$d)
        d_inv <- ifelse(SVD$d > tol, 1 / SVD$d, 0)
        (SVD$v %*% diag(d_inv, length(d_inv)) %*% t(SVD$u)) %*% as.matrix(XtY)
    })
}

# stage-wise estimation: first stage on untreated observations, second stage on residuals
dc_did2s_fit <- function(data, yname, first_stage, second_stage, treatment, weights = NULL) {
    untreat <- data[data[[treatment]] == 0, ]
    w1 <- if (is.null(weights)) NULL else untreat[[weights]]
    f1 <- fixest::xpd(~ 0 + ..fs, lhs = yname, ..fs = first_stage)
    fs <- fixest::feols(f1, data = untreat, weights = w1, fixef.keep_names = TRUE,
                        fixef.rm = "none", warn = FALSE, notes = FALSE)
    first_u <- data[[yname]] - stats::predict(fs, newdata = data)
    data[[yname]] <- first_u
    first_u[data[[treatment]] == 1] <- 0
    w2 <- if (is.null(weights)) NULL else data[[weights]]
    f2 <- fixest::xpd(~ 0 + ..ss, lhs = yname, ..ss = second_stage)
    ss <- fixest::feols(f2, data = data, weights = w2, warn = FALSE, notes = FALSE)
    list(first_stage = fs, second_stage = ss, first_u = first_u)
}

# returns a fixest object whose covariance is the corrected, cluster-robust one
dc_did2s <- function(data, yname, first_stage, second_stage, treatment, cluster_var, weights = NULL) {
    est <- dc_did2s_fit(data, yname, first_stage, second_stage, treatment, weights)
    removed_rows <- est$second_stage$obs_selection$obsRemoved
    if (!is.null(removed_rows)) data <- data[removed_rows, ]
    wv <- if (is.null(weights)) rep.int(1L, nrow(data)) else sqrt(data[[weights]])

    first_u <- est$first_u
    if (!is.null(removed_rows)) first_u <- first_u[removed_rows]
    x1 <- fixest::sparse_model_matrix(est$first_stage, data = data, type = c("rhs", "fixef"))
    second_u <- stats::residuals(est$second_stage)
    x2 <- fixest::sparse_model_matrix(est$second_stage, type = c("rhs", "fixef"))

    first_u <- wv * first_u; x1 <- wv * x1
    second_u <- wv * second_u; x2 <- wv * x2

    x10 <- x1
    treated_rows <- which(data[[treatment]] == 1L) - 1
    x10@x[x10@i %in% treated_rows] <- 0

    x2tx2_inv <- est$second_stage$cov.iid / est$second_stage$sigma2
    IF_ss <- x2tx2_inv %*% Matrix::t(x2 * second_u)
    gamma_hat <- dc_robust_solve(x10, Matrix::crossprod(x1, x2))
    IF_fs <- x2tx2_inv %*% Matrix::t(gamma_hat) %*% Matrix::t(x10 * first_u)
    IF <- IF_fs - IF_ss

    cl <- data[[cluster_var]]
    cov <- Reduce("+", lapply(split(seq_along(cl), cl), function(i)
        Matrix::tcrossprod(Matrix::rowSums(IF[, i, drop = FALSE]))))
    cov <- as.matrix(cov)
    rownames(cov) <- colnames(cov) <- names(est$second_stage$coefficients)
    vl <- list(cov); names(vl) <- sprintf("Corrected Clustered (%s)", cluster_var)
    base::suppressWarnings(summary(est$second_stage, vcov = vl))
}
