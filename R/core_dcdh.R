# Dynamic difference-in-differences estimators for a treatment that can switch on and off
# (de Chaisemartin and D'Haultfoeuille, 2024, Review of Economics and Statistics). Base R only.
# Binary treatment, no covariates. Checked against the R package DIDmultiplegtDYN and Stata did_multiplegt_dyn.
#
# Effect l (l = 1, 2, ...): groups are compared from the period before their first treatment change, F_g - 1, to
# F_g - 1 + l - 1 + 1 = F_g + l - 1, with groups that have the same baseline treatment and have not changed it by then.
# Placebo l: the same comparison l periods BEFORE F_g - 1.

# y, g, t, d: vectors of equal length. Returns the panel as matrices (groups x periods).
dc_dcdh_panel <- function(y, g, t, d) {
    ok <- is.finite(y) & !is.na(g) & !is.na(t) & is.finite(d)
    if (sum(ok) < 10) stop("Too few complete observations.", call. = FALSE)
    if (any(!(d[ok] %in% c(0, 1)))) stop("The treatment must be binary (0/1).", call. = FALSE)
    tt <- sort(unique(t[ok]))
    gg <- sort(unique(g[ok]))
    if (anyDuplicated(cbind(g[ok], t[ok]))) stop("Each group may appear only once per period.", call. = FALSE)
    gi <- match(g[ok], gg); ti <- match(t[ok], tt)
    Y <- matrix(NA_real_, length(gg), length(tt)); D <- Y
    Y[cbind(gi, ti)] <- y[ok]; D[cbind(gi, ti)] <- d[ok]
    list(Y = Y, D = D, groups = gg, times = tt)
}

# Switching structure of each group: baseline treatment, first period with a different treatment, last period observed.
dc_dcdh_design <- function(D) {
    G <- nrow(D); Tn <- ncol(D)
    first <- apply(D, 1, function(z) which(!is.na(z))[1])
    last <- apply(D, 1, function(z) { w <- which(!is.na(z)); w[length(w)] })
    d1 <- D[cbind(seq_len(G), first)]
    Fg <- vapply(seq_len(G), function(i) {
        w <- which(!is.na(D[i, ]) & D[i, ] != d1[i] & seq_len(Tn) > first[i])
        if (length(w)) w[1] else Inf
    }, numeric(1))
    list(d1 = d1, F = Fg, last = last, first = first, sgn = ifelse(is.finite(Fg), ifelse(d1 == 0, 1, -1), NA_real_))
}

# One effect (kind = "effect") or placebo (kind = "placebo") at lag l. Returns the estimate, the number of switchers
# and the influence-function terms u (one per group) from which the variance is computed.
# Switchers with the same baseline treatment and the same first change form a window: they share the periods compared
# and the control groups, and their outcome changes are centred on the window mean (the means are estimated;
# the factor sqrt(n/(n-1)) corrects for it, as in de Chaisemartin and D'Haultfoeuille's variance estimator).
dc_dcdh_one <- function(Y, ds, l, kind = c("effect", "placebo")) {
    kind <- match.arg(kind)
    G <- nrow(Y); Tn <- ncol(Y)
    wins <- list()                                           # one entry per (baseline, F) window
    for (f in sort(unique(ds$F[is.finite(ds$F)]))) for (b in unique(ds$d1[ds$F == f])) {
        sw <- which(ds$F == f & ds$d1 == b)
        t0 <- f - 1
        t1 <- if (kind == "effect") f - 1 + l else f - 1 - l
        if (t0 < 1 || t1 < 1 || t1 > Tn) next
        sw <- sw[!is.na(Y[sw, t0]) & !is.na(Y[sw, t1])]
        if (kind == "effect") sw <- sw[t1 <= ds$last[sw]] else sw <- sw[t1 >= ds$first[sw]]
        cs <- which(ds$d1 == b & ds$F > f - 1 + l & !is.na(Y[, t0]) & !is.na(Y[, t1]))   # not changed by F - 1 + l, also for the placebo
        if (kind == "placebo") {                         # a placebo uses the cells of the matching effect: it needs the outcome at F - 1 + l too
            if (f - 1 + l > Tn) next
            sw <- sw[!is.na(Y[sw, f - 1 + l])]; cs <- cs[!is.na(Y[cs, f - 1 + l])]
        }
        if (length(sw) == 0 || length(cs) == 0) next
        wins[[length(wins) + 1]] <- list(sw = sw, cs = cs, t0 = t0, t1 = t1, s = if (b == 0) 1 else -1)
    }
    n <- sum(vapply(wins, function(w) length(w$sw), numeric(1)))
    u <- numeric(G)
    if (n == 0) return(list(est = NA_real_, n = 0L, u = u))
    est <- 0
    for (w in wins) {
        dS <- Y[w$sw, w$t1] - Y[w$sw, w$t0]
        dC <- Y[w$cs, w$t1] - Y[w$cs, w$t0]
        est <- est + w$s * (sum(dS) - length(w$sw) * mean(dC))
        ns <- length(w$sw); nc <- length(w$cs)
        eS <- if (ns >= 2) mean(dS) else 0; fS <- if (ns >= 2) sqrt(ns / (ns - 1)) else 1
        eC <- if (nc >= 2) mean(dC) else 0; fC <- if (nc >= 2) sqrt(nc / (nc - 1)) else 1
        u[w$sw] <- u[w$sw] + w$s * fS * (dS - eS) / n
        u[w$cs] <- u[w$cs] - w$s * (ns / nc) * fC * (dC - eC) / n
    }
    list(est = est / n, n = as.integer(n), u = u)
}

dc_dcdh <- function(y, g, t, d, effects = 3, placebo = 0, level = 95) {
    pan <- dc_dcdh_panel(y, g, t, d)
    ds <- dc_dcdh_design(pan$D)
    if (!any(is.finite(ds$F))) stop("No group changes its treatment, so no effect can be estimated.", call. = FALSE)
    rows <- list(); U <- NULL
    for (l in seq_len(effects)) rows[[length(rows) + 1]] <- c(list(label = sprintf("Effect_%d", l), l = l, kind = "effect"), dc_dcdh_one(pan$Y, ds, l, "effect"))
    for (l in seq_len(placebo)) rows[[length(rows) + 1]] <- c(list(label = sprintf("Placebo_%d", l), l = l, kind = "placebo"), dc_dcdh_one(pan$Y, ds, l, "placebo"))
    U <- do.call(cbind, lapply(rows, `[[`, "u"))
    V <- crossprod(U)
    est <- vapply(rows, `[[`, numeric(1), "est")
    se <- sqrt(diag(V))
    crit <- stats::qnorm(1 - (1 - level / 100) / 2)
    tab <- data.frame(label = vapply(rows, `[[`, "", "label"), kind = vapply(rows, `[[`, "", "kind"), l = vapply(rows, `[[`, numeric(1), "l"),
                      estimate = est, se = se, lo = est - crit * se, hi = est + crit * se,
                      n_switchers = vapply(rows, `[[`, integer(1), "n"), stringsAsFactors = FALSE)
    tab$se[tab$n_switchers == 0] <- NA
    list(table = tab, vcov = V, G = nrow(pan$Y), n_periods = ncol(pan$Y), n_switchers = sum(is.finite(ds$F)), design = ds, level = level)
}

# Joint test that a set of coefficients is zero.
dc_dcdh_joint <- function(res, kind) {
    w <- which(res$table$kind == kind & is.finite(res$table$estimate) & is.finite(res$table$se))
    if (length(w) == 0) return(NULL)
    b <- res$table$estimate[w]; V <- res$vcov[w, w, drop = FALSE]
    chi <- tryCatch(as.numeric(t(b) %*% solve(V, b)), error = function(e) NA_real_)
    list(chi2 = chi, df = length(w), p = if (is.finite(chi)) stats::pchisq(chi, length(w), lower.tail = FALSE) else NA_real_)
}
