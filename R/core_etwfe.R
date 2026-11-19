# Core layer: Wooldridge's extended two-way fixed effects (ETWFE, Wooldridge 2021, 2023), written without the etwfe package: the formulas of that
# package use i(), which jamovi's formula sandbox blocks. The regression, the covariate treatment and the marginal effects follow the R package
# etwfe (McDermott, 2023), which the tests compare with.
#
# Regression: one dummy for every treated cohort-by-period cell (g, t) with t >= g (control group "not yet treated") or t != g - 1 ("never treated",
# which adds the pre-treatment leads), cell dummies times the covariates demeaned within cohort, covariates, cohort-by-covariate and period-by-
# covariate interactions, and fixed effects for cohort and period (linear model) or cohort and period dummies (Poisson and logit, no fixed effects).
# Effects: the average, over the treated observations of a set, of the change of the prediction when the treatment columns are switched off,
# with delta-method standard errors from the (clustered) covariance matrix of the coefficients.

dc_etwfe_fit <- function(df, covs = character(0), control = c("notyettreated", "nevertreated"), family = c("gaussian", "poisson", "logit")) {
    control <- match.arg(control); family <- match.arg(family)
    y <- df$dc_y; g <- df$dc_g; tt <- df$dc_t
    ug <- sort(unique(g)); ut <- sort(unique(tt)); tref <- min(ut)
    has_never <- any(g == 0)
    if (control == "nevertreated" && !has_never) stop("No never-treated units: choose 'Not yet treated' as control group.", call. = FALSE)
    gref <- if (has_never) 0 else max(ug)
    if (family == "poisson" && any(y < 0)) stop("The Poisson model needs a non-negative outcome.", call. = FALSE)
    if (family == "logit" && !all(y %in% c(0, 1))) stop("The logit model needs an outcome that takes only the values 0 and 1.", call. = FALSE)
    cohorts <- setdiff(ug, gref)
    if (!length(cohorts)) stop("There are no treated cohorts.", call. = FALSE)
    # cohort means of the covariates (taken over all observations, as etwfe does before it drops rows)
    X <- if (length(covs)) as.matrix(df[, covs, drop = FALSE]) else NULL
    Xdm <- if (!is.null(X)) apply(X, 2, function(v) v - stats::ave(v, g)) else NULL
    if (!is.null(Xdm)) Xdm <- matrix(Xdm, nrow = nrow(df), dimnames = list(NULL, covs))
    keep <- rep(TRUE, nrow(df))
    if (control == "notyettreated" && !has_never) keep <- tt < gref                                   # the last cohort is the control until it is treated
    # treated cells
    # the first period is left out of the cells with not-yet-treated controls only (etwfe: ref2 = tref); with never-treated controls every period but g - 1 is a cell
    cells <- do.call(rbind, lapply(cohorts, function(G) data.frame(g = G, t = if (control == "notyettreated") ut[ut != tref & ut >= G] else ut[ut != G - 1])))
    cells <- cells[cells$t %in% unique(tt[keep]), , drop = FALSE]
    if (control == "notyettreated" && !has_never) cells <- cells[cells$g < gref, , drop = FALSE]
    if (!nrow(cells)) stop("There are no treated cohort-by-period cells to estimate.", call. = FALSE)
    nm_cell <- sprintf("dc_c_%d", seq_len(nrow(cells)))
    M <- matrix(0, nrow(df), nrow(cells), dimnames = list(NULL, nm_cell))
    for (k in seq_len(nrow(cells))) M[, k] <- as.numeric(g == cells$g[k] & tt == cells$t[k])
    treat_cols <- nm_cell
    if (!is.null(Xdm)) for (j in seq_along(covs)) {
        Mj <- M[, nm_cell, drop = FALSE] * Xdm[, j]; colnames(Mj) <- sprintf("dc_cx%d_%d", j, seq_len(nrow(cells)))
        M <- cbind(M, Mj); treat_cols <- c(treat_cols, colnames(Mj))
    }
    if (!is.null(X)) {
        M <- cbind(M, X); colnames(M)[ncol(M) - ncol(X) + seq_len(ncol(X))] <- sprintf("dc_x_%d", seq_along(covs))
        for (j in seq_along(covs)) {
            for (G in setdiff(ug, gref)) { v <- as.numeric(g == G) * X[, j]; M <- cbind(M, v); colnames(M)[ncol(M)] <- sprintf("dc_gx%d_%s", j, gsub("[^0-9A-Za-z]", "_", format(G))) }
            for (T in setdiff(ut, tref)) { v <- as.numeric(tt == T) * X[, j]; M <- cbind(M, v); colnames(M)[ncol(M)] <- sprintf("dc_tx%d_%s", j, gsub("[^0-9A-Za-z]", "_", format(T))) }
        }
    }
    dat <- data.frame(M, check.names = FALSE); dat$dc_y <- y; dat$dc_cl <- df$dc_cl
    if (family == "gaussian") {
        dat$dc_gf <- factor(g); dat$dc_tf <- factor(tt)
    } else {
        for (G in setdiff(ug, gref)) { dat[[sprintf("dc_g_%s", gsub("[^0-9A-Za-z]", "_", format(G)))]] <- as.numeric(g == G) }
        for (T in setdiff(ut, tref)) { dat[[sprintf("dc_t_%s", gsub("[^0-9A-Za-z]", "_", format(T)))]] <- as.numeric(tt == T) }
    }
    dat <- dat[keep, , drop = FALSE]
    # drop columns that are constant zero (cells without observations) before fitting
    rhs <- setdiff(names(dat), c("dc_y", "dc_cl", "dc_gf", "dc_tf"))
    rhs <- rhs[vapply(rhs, function(v) any(dat[[v]] != 0), TRUE)]
    fml <- if (family == "gaussian") stats::as.formula(paste("dc_y ~", paste(rhs, collapse = " + "), "| dc_gf + dc_tf"))
           else stats::as.formula(paste("dc_y ~", paste(rhs, collapse = " + ")))
    fit <- if (family == "gaussian") fixest::feols(fml, data = dat, vcov = ~dc_cl, notes = FALSE, nthreads = 1)
           else fixest::feglm(fml, data = dat, family = if (family == "poisson") "poisson" else "binomial", vcov = ~dc_cl, notes = FALSE, nthreads = 1)
    b <- stats::coef(fit); V <- stats::vcov(fit)
    b <- b[!is.na(b)]; V <- V[names(b), names(b), drop = FALSE]
    list(fit = fit, b = b, V = V, dat = dat, keep = keep, cells = cells, treat_cols = intersect(treat_cols, names(b)), family = family, control = control,
         gref = gref, tref = tref, g = g[keep], t = tt[keep], has_never = has_never, covs = covs, nm_cell = nm_cell)
}

# Effect of each treated observation and its gradient with respect to the coefficients.
dc_etwfe_obs <- function(m) {
    nb <- names(m$b); dat <- m$dat
    Xf <- as.matrix(dat[, intersect(nb, names(dat)), drop = FALSE])
    if ("(Intercept)" %in% nb) Xf <- cbind(`(Intercept)` = 1, Xf)
    Xf <- Xf[, nb, drop = FALSE]
    Z <- Xf; Z[, setdiff(nb, m$treat_cols)] <- 0                                     # only the treatment columns
    tau <- as.numeric(Z %*% m$b)
    if (m$family == "gaussian") return(list(effect = tau, grad = Z, active = rowSums(Z != 0) > 0 | rowSums(as.matrix(dat[, m$nm_cell, drop = FALSE])) > 0))
    Gf <- if (m$family == "poisson") list(f = exp, d = exp) else list(f = stats::plogis, d = function(e) stats::plogis(e) * (1 - stats::plogis(e)))
    eta <- as.numeric(Xf %*% m$b); X0 <- Xf - Z; eta0 <- eta - tau
    list(effect = Gf$f(eta) - Gf$f(eta0), grad = Gf$d(eta) * Xf - Gf$d(eta0) * X0, active = rowSums(as.matrix(dat[, m$nm_cell, drop = FALSE])) > 0)
}

# Aggregation of the effects over the treated observations (post-treatment ones, or all cells for the event study of the never-treated design).
# type: "simple", "group", "calendar" or "event". Returns a data frame (key, est, se, n).
dc_etwfe_agg <- function(m, type = c("simple", "group", "calendar", "event")) {
    type <- match.arg(type)
    o <- dc_etwfe_obs(m)
    g <- m$g; tt <- m$t
    sel <- o$active & g != m$gref
    if (type != "event") sel <- sel & tt >= g
    key <- switch(type, simple = rep(1, length(g)), group = g, calendar = tt, event = tt - g)
    idx <- which(sel)
    out <- lapply(sort(unique(key[idx])), function(k) {
        i <- idx[key[idx] == k]
        gr <- colMeans(o$grad[i, , drop = FALSE]); est <- mean(o$effect[i])
        data.frame(key = k, est = est, se = sqrt(max(drop(t(gr) %*% m$V %*% gr), 0)), n = length(i))
    })
    do.call(rbind, out)
}

# Estimator adapter in the common format of core_estimators.R.
dc_etwfe <- function(d, min_e, max_e, control = c("notyettreated", "nevertreated"), family = c("gaussian", "poisson", "logit")) {
    control <- match.arg(control); family <- match.arg(family)
    method <- paste0("Extended TWFE (Wooldridge)", if (family != "gaussian") paste0(", ", family) else "")
    tryCatch({
        d <- dc_balanced_subset_etwfe(d)
        cx <- d$info$covs
        covs <- if (is.null(cx) || !nrow(cx)) character(0) else cx$code
        m <- dc_etwfe_fit(d$data, covs, control, family)
        simple <- dc_etwfe_agg(m, "simple"); ev <- dc_etwfe_agg(m, "event")
        grp <- dc_etwfe_agg(m, "group"); cal <- dc_etwfe_agg(m, "calendar")
        dyn <- dc_dyn(ev$key, ev$est, ev$se, min_e, max_e)
        dc_out(method, att = simple$est[1], se = simple$se[1], dynamic = dyn,
               control = if (control == "nevertreated") "Never-treated units" else "Not-yet-treated and never-treated units",
               extra = list(group = data.frame(g = grp$key, est = grp$est, se = grp$se), calendar = data.frame(t = cal$key, est = cal$est, se = cal$se)),
               n_obs = nrow(m$dat), n_cl = length(unique(m$dat$dc_cl)))
    }, error = function(e) dc_fail(method, e))
}

# ETWFE does not need a balanced panel; the data are used as they are.
dc_balanced_subset_etwfe <- function(d) d
