# Core layer: diagnostics and sensitivity analysis. No jamovi dependencies.

dc_wald <- function(b, V) {
    k <- length(b)
    stat <- drop(t(b) %*% solve(V, b))
    list(stat = stat, df = k, p = stats::pchisq(stat, k, lower.tail = FALSE))
}

# ---- Goodman-Bacon decomposition of the static TWFE coefficient -------------------------
# Computed from cohort-level means, which is exact for balanced panels and takes milliseconds
# instead of the minutes bacondecomp::bacon() needs on large panels (it fits one regression with
# unit dummies per 2x2 comparison). The formulas are those of Goodman-Bacon (2021) as implemented
# in bacondecomp (group sizes are unit-periods inside each comparison window); dev/test_core.R
# checks that weights and estimates are identical.
# Returns the 2x2 comparisons (`detail`: treated, untreated, estimate, weight, type), a summary
# by comparison type and the TWFE coefficient implied by the decomposition (`twfe`).
dc_bacon_2x2 <- function(df) {
    never <- 99999
    g <- ifelse(df$dc_g == 0, never, df$dc_g)
    times <- sort(unique(df$dc_t))
    groups <- sort(unique(g))
    n_k <- tapply(!duplicated(df$dc_id), g, sum)                  # units per cohort
    ybar <- tapply(df$dc_y, list(factor(g, levels = groups), factor(df$dc_t, levels = times)), mean)
    out <- list()
    for (k in groups[groups != never & groups != min(times)]) {
        for (l in groups[groups != k]) {
            w <- if (l == never) times else if (k < l) times[times < l] else times[times >= l]
            treated_in <- function(gr) mean(w >= gr)                 # share of the window treated
            rows <- function(gr) n_k[[as.character(gr)]] * length(w)   # unit-periods of the group inside the window
            if (l == never) {
                nk <- rows(k); nu <- rows(never)
                nku <- nk / (nk + nu); Dk <- treated_in(k)
                weight <- (nk + nu)^2 * (nku * (1 - nku) * Dk * (1 - Dk))
                type <- "Treated vs Untreated"
            } else if (k < l) {
                nk <- rows(k); nl <- rows(l)
                nkl <- nk / (nk + nl); Dk <- treated_in(k); Dl <- treated_in(l)
                weight <- ((nk + nl) * (1 - Dl))^2 *
                    (nkl * (1 - nkl) * (Dk - Dl) / (1 - Dl) * (1 - Dk) / (1 - Dl))
                type <- "Earlier vs Later Treated"
            } else {
                nk <- rows(l); nl <- rows(k)
                nkl <- nk / (nk + nl); Dk <- treated_in(l); Dl <- treated_in(k)
                weight <- ((nk + nl) * Dk)^2 * (nkl * (1 - nkl) * (Dl / Dk) * (Dk - Dl) / Dk)
                type <- "Later vs Earlier Treated"
            }
            pre <- w[w < k]; post <- w[w >= k]
            m <- function(gr, tt) mean(ybar[as.character(gr), as.character(tt)])
            est <- (m(k, post) - m(k, pre)) - (m(l, post) - m(l, pre))
            out[[length(out) + 1]] <- data.frame(treated = k, untreated = l, estimate = est,
                                                 weight = weight, type = type, stringsAsFactors = FALSE)
        }
    }
    res <- do.call(rbind, out)
    res$weight <- res$weight / sum(res$weight)
    res
}

dc_bacon <- function(d) {
    tryCatch({
        if (!d$info$balanced)
            stop("The Goodman-Bacon decomposition requires a balanced panel.")
        x <- dc_bacon_2x2(d$data)
        types <- unique(x$type)
        summ <- do.call(rbind, lapply(types, function(ty) {
            s <- x[x$type == ty, ]
            data.frame(type = ty, n = nrow(s), weight = sum(s$weight),
                       est = sum(s$weight * s$estimate) / sum(s$weight),
                       stringsAsFactors = FALSE)
        }))
        list(ok = TRUE, message = "", detail = x, summary = summ,
             twfe = sum(x$weight * x$estimate))
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
}

# ---- TWFE weights on treated observations (de Chaisemartin & D'Haultfoeuille 2020) -----
# beta_TWFE = sum over treated cells of w_it * tau_it, with w_it proportional to the
# residual of the treatment indicator after removing unit and period effects.
dc_twfe_weights <- function(d) {
    tryCatch({
        df <- d$data
        m <- fixest::feols(dc_d ~ 1 | dc_id + dc_t, data = df, notes = FALSE, nthreads = 1)
        e <- as.numeric(stats::resid(m))
        tr <- df$dc_d == 1
        w <- e[tr] / sum(e[tr])
        g <- df$dc_g[tr]
        by_g <- do.call(rbind, lapply(sort(unique(g)), function(k) {
            wk <- w[g == k]
            data.frame(cohort = k, n = length(wk), weight = sum(wk),
                       neg_n = sum(wk < 0), neg_weight = sum(wk[wk < 0]),
                       stringsAsFactors = FALSE)
        }))
        list(ok = TRUE, message = "", w = w,
             n_treated = length(w), pos_n = sum(w > 0), neg_n = sum(w < 0),
             pos_w = sum(w[w > 0]), neg_w = sum(w[w < 0]), by_cohort = by_g)
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
}

# ---- Joint Wald tests of the pre-treatment (e < -1) coefficients ------------------------
# Returns a data frame with columns method, stat, df, p, note.
dc_pretests <- function(d, min_e, control = c("nevertreated", "notyettreated")) {
    control <- match.arg(control)
    one <- function(method, f) {
        tryCatch(f(), error = function(e)
            data.frame(method = method, stat = NA_real_, df = NA_real_, p = NA_real_,
                       note = conditionMessage(e), stringsAsFactors = FALSE))
    }
    row <- function(method, w, note = "")
        data.frame(method = method, stat = w$stat, df = w$df, p = w$p, note = note,
                   stringsAsFactors = FALSE)
    df <- d$data
    out <- list()

    out$twfe <- one("TWFE", function() {
        me <- fixest::feols(dc_y ~ i(dc_rel, ref = c(-1, -1000)) | dc_id + dc_t,
                            data = df, vcov = ~dc_id, notes = FALSE, nthreads = 1)
        b <- stats::coef(me); e <- as.numeric(sub(".*::", "", names(b)))
        k <- which(e < -1 & e >= min_e)
        if (!length(k)) stop("No pre-treatment periods in the event window.")
        row("TWFE", dc_wald(b[k], stats::vcov(me)[k, k, drop = FALSE]))
    })

    out$sa <- one("Sun-Abraham", function() {
        core <- dc_sa_core(d)
        k <- which(core$e < -1 & core$e >= min_e)
        if (!length(k)) stop("No pre-treatment periods in the event window.")
        row("Sun-Abraham", dc_wald(core$est[k], core$V[k, k, drop = FALSE]))
    })

    out$cs <- one("Callaway-Sant'Anna", function() {
        d <- dc_balanced_subset(d); df <- d$data
        if (control == "nevertreated" && d$info$n_never == 0)
            stop("No never-treated units: choose 'Not yet treated' as control group.")
        fit <- did::att_gt(yname = "dc_y", tname = "dc_t", idname = "dc_id", gname = "dc_g",
                           data = df, control_group = control, base_period = "universal",
                           bstrap = FALSE, cband = FALSE, print_details = FALSE)
        if (is.null(fit$Wpval) || is.na(fit$Wpval))
            stop("The pre-test is not available (no pre-treatment group-time effects).")
        npre <- sum(fit$group > fit$t + 1)
        row("Callaway-Sant'Anna", list(stat = fit$W, df = npre, p = fit$Wpval),
            "Tests all pre-treatment ATT(g,t) jointly")
    })
    do.call(rbind, out)
}

# ---- HonestDiD sensitivity (Rambachan & Roth 2023), relative magnitudes -----------------
# Prepares event-study estimates and their covariance (pre e < -1, post e >= 0; the
# reference period e = -1 is dropped). Pre and post periods must be contiguous.
dc_honest_prepare <- function(d, base = c("cs", "sa"), min_e, max_e,
                              cs_control = "nevertreated") {
    base <- match.arg(base)
    if (base == "sa") {
        core <- dc_sa_core(d)
        e <- core$est; ee <- core$e; V <- core$V
        keep <- ee != -1 & ee >= min_e & ee <= max_e
        beta <- core$est[keep]; sigma <- core$V[keep, keep, drop = FALSE]; e <- ee[keep]
    } else {
        d <- dc_balanced_subset(d)
        if (cs_control == "nevertreated" && d$info$n_never == 0)
            stop("No never-treated units: choose 'Not yet treated' as control group.")
        fit <- did::att_gt(yname = "dc_y", tname = "dc_t", idname = "dc_id", gname = "dc_g",
                           data = d$data, control_group = cs_control, base_period = "universal",
                           bstrap = FALSE, cband = FALSE, print_details = FALSE)
        es <- did::aggte(fit, type = "dynamic", min_e = min_e, max_e = max_e,
                         bstrap = FALSE, cband = FALSE)
        IF <- es$inf.function$dynamic.inf.func.e
        n <- nrow(IF)
        V <- crossprod(IF) / n / n
        keep <- es$egt != -1
        beta <- es$att.egt[keep]; sigma <- V[keep, keep, drop = FALSE]; e <- es$egt[keep]
    }
    pre <- e[e < -1]; post <- e[e >= 0]
    if (length(pre) < 1)
        stop("HonestDiD needs at least one pre-treatment period besides e = -1 (widen the event window).")
    if (length(post) < 1)
        stop("HonestDiD needs at least one post-treatment period.")
    contiguous <- function(x, first, last) isTRUE(all.equal(as.numeric(x), as.numeric(seq(first, last))))
    if (!contiguous(pre, min(pre), -2) || !contiguous(post, 0, max(post)))
        stop("HonestDiD needs consecutive pre- and post-treatment event times; some are missing in the window.")
    list(beta = beta, sigma = sigma, npre = length(pre), npost = length(post), e = e)
}

# One point of the sensitivity analysis. mbar = NA gives the original (unadjusted) interval.
dc_honest_one <- function(prep, mbar, target = c("first", "average"), alpha = 0.05) {
    target <- match.arg(target)
    l <- if (target == "first") hd_basisVector(1, prep$npost) else rep(1 / prep$npost, prep$npost)
    res <- suppressWarnings(
        if (is.na(mbar))
            hd_constructOriginalCS(prep$beta, prep$sigma, prep$npre, prep$npost,
                                           l_vec = l, alpha = alpha)
        else
            hd_createSensitivityResults_relativeMagnitudes(
                prep$beta, prep$sigma, numPrePeriods = prep$npre, numPostPeriods = prep$npost,
                Mbarvec = mbar, l_vec = l, alpha = alpha))
    c(lb = as.numeric(res$lb[1]), ub = as.numeric(res$ub[1]))
}
