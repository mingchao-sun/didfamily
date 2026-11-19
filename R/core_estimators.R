# Core layer: estimator adapters. Each returns the same structure:
#   list(method, ok, message, att, se, dynamic = data.frame(e, est, se), control)
# `dynamic` always contains the reference period e = -1 (estimate 0, se NA).

# `df`: degrees of freedom of the t distribution used for inference (NA = normal approximation).
# Regression-based estimators (TWFE, Sun-Abraham, stacked) use t(G - 1) with G clusters, like fixest and
# Stata's reghdfe; Callaway-Sant'Anna, imputation and Gardner keep their packages' normal inference.
dc_out <- function(method, att = NA_real_, se = NA_real_, dynamic = NULL,
                   control = "", ok = TRUE, message = "", extra = NULL,
                   df = NA_real_, n_obs = NA_integer_, n_cl = NA_integer_) {
    list(method = method, ok = ok, message = message,
         att = att, se = se, dynamic = dynamic, control = control, extra = extra,
         df = df, n_obs = n_obs, n_cl = n_cl)
}

# critical value, interval and p-value for an estimate under the estimator's own reference distribution
dc_inference <- function(est, se, df = NA_real_, level = 95) {
    a <- 1 - level / 100
    if (is.finite(df)) {
        crit <- stats::qt(1 - a / 2, df); p <- 2 * stats::pt(-abs(est / se), df)
    } else {
        crit <- stats::qnorm(1 - a / 2); p <- 2 * stats::pnorm(-abs(est / se))
    }
    list(crit = crit, lo = est - crit * se, hi = est + crit * se, p = p)
}

# sampling weights as a vector, or NULL when the analysis is unweighted
dc_sw <- function(d) if (isTRUE(d$info$weighted)) d$data$dc_sw else NULL

dc_fail <- function(method, e) {
    dc_out(method, ok = FALSE, message = conditionMessage(e))
}

# add the e = -1 reference row, clip to window, sort
dc_dyn <- function(e, est, se, min_e, max_e) {
    dyn <- data.frame(e = as.numeric(e), est = as.numeric(est), se = as.numeric(se))
    dyn <- dyn[dyn$e != -1, , drop = FALSE]
    dyn <- rbind(dyn, data.frame(e = -1, est = 0, se = NA_real_))
    dyn <- dyn[dyn$e >= min_e & dyn$e <= max_e, , drop = FALSE]
    dyn[order(dyn$e), , drop = FALSE]
}

# Covariate terms as a language object (NULL when there are none). Formulas are built from
# code templates plus the internal column names dc_x1, dc_x2, ... (never from user text).
# `fe = TRUE` (default) keeps only covariates that vary within units: the others are
# collinear with unit fixed effects. Callaway-Sant'Anna and DDD pass fe = FALSE.
dc_cov_expr <- function(d, fe = TRUE) {
    cm <- d$info$covs
    codes <- if (is.null(cm)) character() else if (fe) cm$code[cm$varying] else cm$code
    if (length(codes)) str2lang(paste(codes, collapse = " + ")) else NULL
}

dc_twfe <- function(d, min_e, max_e) {
    method <- "TWFE"
    tryCatch({
        df <- d$data
        cx <- dc_cov_expr(d)
        f_static <- if (is.null(cx)) dc_y ~ dc_d | dc_id + dc_t
                    else eval(bquote(dc_y ~ dc_d + .(cx) | dc_id + dc_t))
        f_event <- if (is.null(cx)) dc_y ~ i(dc_rel, ref = c(-1, -1000)) | dc_id + dc_t
                   else eval(bquote(dc_y ~ i(dc_rel, ref = c(-1, -1000)) + .(cx) | dc_id + dc_t))
        m <- fixest::feols(f_static, data = df, weights = dc_sw(d), vcov = ~dc_cl, notes = FALSE, nthreads = 1)
        ct <- fixest::coeftable(m)
        me <- fixest::feols(f_event, data = df, weights = dc_sw(d), vcov = ~dc_cl, notes = FALSE, nthreads = 1)
        ce <- fixest::coeftable(me)
        ce <- ce[grepl("^dc_rel::", rownames(ce)), , drop = FALSE]
        e <- as.numeric(sub(".*::", "", rownames(ce)))
        dc_out(method, att = ct["dc_d", 1], se = ct["dc_d", 2],
               dynamic = dc_dyn(e, ce[, 1], ce[, 2], min_e, max_e),
               control = "Never + not-yet treated (implicit weights)",
               df = fixest::degrees_freedom(m, "t"), n_obs = stats::nobs(m), n_cl = length(unique(df$dc_cl)))
    }, error = function(e) dc_fail(method, e))
}

# Sun & Abraham (2021) interaction-weighted estimator.
# Implemented directly instead of fixest::sunab(): jamovi replaces stats::as.formula
# with a validator that rejects sunab(), while plain i() interactions are fine.
# Cohort x event-time effects are estimated against never-treated units, then
# averaged with weights proportional to cohort sizes (treated as fixed, delta method).
# Returns event-time estimates `est` with full covariance `V` (all e except -1),
# plus the overall ATT over post-treatment cells.
dc_sa_core <- function(d) {
    df <- d$data
    if (d$info$n_never == 0)
        stop("Sun-Abraham in this module needs a never-treated group.")
    df$dc_gf <- factor(df$dc_g)   # factor so i() builds one coefficient per cohort
    cx <- dc_cov_expr(d)
    f <- if (is.null(cx)) dc_y ~ i(dc_rel, dc_gf, ref = -1, ref2 = "0") | dc_id + dc_t
         else eval(bquote(dc_y ~ i(dc_rel, dc_gf, ref = -1, ref2 = "0") + .(cx) | dc_id + dc_t))
    m <- fixest::feols(f, data = df, weights = dc_sw(d), vcov = ~dc_cl, notes = FALSE, nthreads = 1)
    b <- stats::coef(m)
    V <- stats::vcov(m)
    mt <- regmatches(names(b), regexec("^dc_rel::(-?[0-9]+):dc_gf::([0-9]+)$", names(b)))
    ok <- lengths(mt) == 3
    if (!any(ok))
        stop("Could not identify cohort-by-event-time coefficients.")
    cells <- data.frame(idx = which(ok),
                        e = as.numeric(vapply(mt[ok], `[`, "", 2)),
                        g = as.numeric(vapply(mt[ok], `[`, "", 3)))
    # cell sizes for the aggregation weights: unit counts, or sums of the unit weights when weighted
    n_cell <- if (isTRUE(d$info$weighted)) tapply(df$dc_sw, list(factor(df$dc_g), factor(df$dc_rel)), sum)
              else table(df$dc_g, df$dc_rel)
    cells$n <- n_cell[cbind(as.character(cells$g), as.character(cells$e))]
    bc <- b[cells$idx]
    Vc <- V[cells$idx, cells$idx, drop = FALSE]

    es <- sort(unique(cells$e))
    Wm <- t(vapply(es, function(x) {
        w <- numeric(nrow(cells)); i <- cells$e == x
        w[i] <- cells$n[i] / sum(cells$n[i]); w
    }, numeric(nrow(cells))))
    wa <- ifelse(cells$e >= 0, cells$n / sum(cells$n[cells$e >= 0]), 0)
    V_e <- Wm %*% Vc %*% t(Wm)
    list(e = es, est = drop(Wm %*% bc), V = V_e,
         att = sum(wa * bc), se_att = sqrt(drop(t(wa) %*% Vc %*% wa)),
         df = fixest::degrees_freedom(m, "t"), n_obs = stats::nobs(m))
}

dc_sa <- function(d, min_e, max_e) {
    method <- "Sun-Abraham"
    tryCatch({
        core <- dc_sa_core(d)
        dc_out(method, att = core$att, se = core$se_att,
               dynamic = dc_dyn(core$e, core$est, sqrt(diag(core$V)), min_e, max_e),
               control = "Never treated", df = core$df, n_obs = core$n_obs,
               n_cl = length(unique(d$data$dc_cl)))
    }, error = function(e) dc_fail(method, e))
}

# `boot`: list(on, reps, seed) for the multiplier bootstrap. A cluster variable above the
# unit level also switches to the bootstrap (did only supports clustering that way).
# `extras = TRUE` adds group-time effects and group / calendar-time aggregations.
dc_cs <- function(d, min_e, max_e, control = c("nevertreated", "notyettreated"),
                  est_method = c("dr", "ipw", "reg"), boot = NULL, extras = FALSE) {
    method <- "Callaway-Sant'Anna"
    control <- match.arg(control)
    est_method <- match.arg(est_method)
    # the bootstrap seeds the generator; put the caller's random-number state back afterwards
    if (isTRUE(d$info$cluster) || isTRUE(boot$on)) {
        had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
        old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
        on.exit({
            if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
        }, add = TRUE)
    }
    tryCatch({
        d <- dc_balanced_subset(d)
        df <- d$data
        if (control == "nevertreated" && d$info$n_never == 0)
            stop("No never-treated units: choose 'Not yet treated' as control group.")
        cx <- dc_cov_expr(d, fe = FALSE)
        xf <- if (is.null(cx)) ~1 else eval(bquote(~ .(cx)))
        clustered <- isTRUE(d$info$cluster)
        use_boot <- clustered || isTRUE(boot$on)
        reps <- if (is.null(boot$reps)) 999L else as.integer(boot$reps)
        seed <- if (is.null(boot$seed)) 1234L else as.integer(boot$seed)
        cv <- if (clustered) "dc_cl" else NULL
        agg <- function(fit, type, ...) {
            if (use_boot) set.seed(seed)
            did::aggte(fit, type = type, bstrap = use_boot, biters = reps, cband = FALSE,
                       clustervars = cv, ...)
        }
        if (use_boot) set.seed(seed)
        fit <- did::att_gt(yname = "dc_y", tname = "dc_t", idname = "dc_id",
                           gname = "dc_g", xformla = xf, data = df, control_group = control,
                           est_method = est_method,
                           weightsname = if (isTRUE(d$info$weighted)) "dc_sw" else NULL,
                           base_period = "universal", bstrap = use_boot, biters = reps,
                           clustervars = cv, cband = FALSE, print_details = FALSE)
        simple <- agg(fit, "simple")
        dyn <- agg(fit, "dynamic", min_e = min_e, max_e = max_e)
        extra <- NULL
        if (extras) {
            grp <- agg(fit, "group"); cal <- agg(fit, "calendar")
            extra <- list(
                gt = data.frame(g = fit$group, t = fit$t, att = fit$att, se = as.numeric(fit$se)),
                group = data.frame(g = grp$egt, est = grp$att.egt, se = grp$se.egt),
                calendar = data.frame(t = cal$egt, est = cal$att.egt, se = cal$se.egt))
        }
        ctl <- if (control == "nevertreated") "Never treated" else "Not yet treated"
        if (use_boot) ctl <- paste0(ctl, "; bootstrap SE (", reps, " draws)")
        if (!is.null(d$info$n_unbalanced_dropped))
            ctl <- paste0(ctl, sprintf("; %d units with missing periods set aside", d$info$n_unbalanced_dropped))
        dc_out(method, att = simple$overall.att, se = simple$overall.se,
               dynamic = dc_dyn(dyn$egt, dyn$att.egt, dyn$se.egt, min_e, max_e),
               control = ctl, extra = extra, n_obs = nrow(df), n_cl = length(unique(df$dc_cl)))
    }, error = function(e) dc_fail(method, e))
}

# Borusyak, Jaravel & Spiess (2024) imputation estimator via didimputation.
# Event-study mode of did_imputation() builds an i() formula that jamovi's formula
# validator rejects, so post-treatment horizons are obtained one at a time with
# aggregation weights (weight 1 on untreated rows keeps the first stage intact).
# Pre-treatment coefficients follow BJS's pre-trend test: lead dummies are added to a
# regression fitted on untreated observations only (never-treated units and treated units
# before treatment), so they test the parallel-trends assumption behind the imputation.
dc_bjs <- function(d, min_e, max_e) {
    method <- "Borusyak-Jaravel-Spiess"
    tryCatch({
        if (isTRUE(d$info$weighted))
            stop("didimputation has no sampling weights, so Borusyak-Jaravel-Spiess is skipped when a weight variable is used.")
        df <- d$data
        cx <- dc_cov_expr(d)
        fs <- if (is.null(cx)) NULL else eval(bquote(~ .(cx) | dc_id + dc_t))
        run <- function(w = NULL) {
            if (!is.null(w)) df$dc_w <- w
            r <- didimputation::did_imputation(
                yname = "dc_y", gname = "dc_g", tname = "dc_t", idname = "dc_id",
                data = df, first_stage = fs, wname = if (is.null(w)) NULL else "dc_w",
                cluster_var = "dc_cl")
            r <- as.data.frame(r)
            c(est = r$estimate[1], se = r$std.error[1])
        }
        all <- run()
        post <- sort(unique(df$dc_rel[df$dc_d == 1 & df$dc_rel >= 0 & df$dc_rel <= max_e]))
        by_e <- t(vapply(post, function(h)
            run(ifelse(df$dc_d == 1, as.numeric(df$dc_rel == h), 1)), c(est = 0, se = 0)))
        est <- by_e[, "est"]; se <- by_e[, "se"]

        # pre-treatment coefficients (leads), from untreated observations only
        leads <- sort(unique(df$dc_rel[df$dc_g > 0 & df$dc_rel < -1 & df$dc_rel >= min_e]))
        if (length(leads)) {
            u <- df[df$dc_d == 0, , drop = FALSE]
            nm <- paste0("dc_ld_", abs(leads))
            for (k in seq_along(leads)) u[[nm[k]]] <- as.numeric(u$dc_g > 0 & u$dc_rel == leads[k])
            rhs <- paste(c(nm, if (!is.null(cx)) deparse(cx)), collapse = " + ")
            f <- eval(bquote(dc_y ~ .(str2lang(rhs)) | dc_id + dc_t))
            m <- fixest::feols(f, data = u, vcov = ~dc_cl, notes = FALSE, nthreads = 1)
            ct <- fixest::coeftable(m)
            ct <- ct[rownames(ct) %in% nm, , drop = FALSE]
            lead_e <- leads[match(rownames(ct), nm)]
            est <- c(est, ct[, 1]); se <- c(se, ct[, 2]); post <- c(post, lead_e)
        }
        dc_out(method, att = all[["est"]], se = all[["se"]],
               dynamic = dc_dyn(post, est, se, min_e, max_e),
               control = "Untreated observations (imputed counterfactual)",
               n_obs = nrow(df), n_cl = length(unique(df$dc_cl)))
    }, error = function(e) dc_fail(method, e))
}

# Gardner (2022) two-stage DiD via did2s. Event-time dummies are built by hand
# (did2s's i() interface is blocked by jamovi's formula validator).
dc_gardner <- function(d, min_e, max_e) {
    method <- "Gardner (did2s)"
    tryCatch({
        df <- d$data
        trt <- df$dc_g > 0
        es <- sort(setdiff(unique(df$dc_rel[trt]), -1))
        nm <- ifelse(es < 0, paste0("dc_ev_m", abs(es)), paste0("dc_ev_p", es))
        for (k in seq_along(es))
            df[[nm[k]]] <- as.numeric(trt & df$dc_rel == es[k])

        cx <- dc_cov_expr(d)
        fs <- if (is.null(cx)) ~ 0 | dc_id + dc_t else eval(bquote(~ .(cx) | dc_id + dc_t))
        wn <- if (isTRUE(d$info$weighted)) "dc_sw" else NULL
        static <- dc_did2s(df, yname = "dc_y", first_stage = fs,
                               second_stage = ~ dc_d, treatment = "dc_d",
                               cluster_var = "dc_cl", weights = wn)
        cs <- fixest::coeftable(static)
        dyn <- dc_did2s(df, yname = "dc_y", first_stage = fs,
                            second_stage = stats::as.formula(paste("~", paste(nm, collapse = " + "))),
                            treatment = "dc_d", cluster_var = "dc_cl", weights = wn)
        cd <- fixest::coeftable(dyn)
        e <- es[match(rownames(cd), nm)]
        dc_out(method, att = cs["dc_d", 1], se = cs["dc_d", 2],
               dynamic = dc_dyn(e, cd[, 1], cd[, 2], min_e, max_e),
               control = "Untreated observations (imputed counterfactual)",
               n_obs = nrow(df), n_cl = length(unique(df$dc_cl)))
    }, error = function(e) dc_fail(method, e))
}

# Stacked DiD (Cengiz et al. 2019): one sub-experiment per treatment cohort, built
# from the cohort plus clean controls (never treated, or not treated by the end of the
# window), stacked with unit-by-experiment and period-by-experiment fixed effects.
dc_stacked <- function(d, min_e, max_e) {
    method <- "Stacked DiD"
    tryCatch({
        df <- d$data
        cohorts <- sort(unique(df$dc_g[df$dc_g > 0]))
        stacks <- lapply(cohorts, function(g) {
            hi <- g + max_e
            w <- df[df$dc_t >= g + min_e & df$dc_t <= hi, , drop = FALSE]
            w <- w[w$dc_g == g | w$dc_g == 0 | w$dc_g > hi, , drop = FALSE]
            if (!any(w$dc_g != g)) return(NULL)      # no clean control for this cohort
            w$dc_trt   <- as.numeric(w$dc_g == g)
            w$dc_erel  <- w$dc_t - g
            w$dc_trtpost <- w$dc_trt * as.numeric(w$dc_erel >= 0)
            w$dc_stack <- g
            w
        })
        stacks <- stacks[!vapply(stacks, is.null, TRUE)]
        if (!length(stacks))
            stop("No cohort has clean controls inside the event window.")
        st <- do.call(rbind, stacks)
        st$dc_ue <- paste(st$dc_id, st$dc_stack)
        st$dc_te <- paste(st$dc_t, st$dc_stack)

        cx <- dc_cov_expr(d)
        f_static <- if (is.null(cx)) dc_y ~ dc_trtpost | dc_ue + dc_te
                    else eval(bquote(dc_y ~ dc_trtpost + .(cx) | dc_ue + dc_te))
        f_event <- if (is.null(cx)) dc_y ~ i(dc_erel, dc_trt, ref = -1) | dc_ue + dc_te
                   else eval(bquote(dc_y ~ i(dc_erel, dc_trt, ref = -1) + .(cx) | dc_ue + dc_te))
        sw <- if (isTRUE(d$info$weighted)) st$dc_sw else NULL
        m <- fixest::feols(f_static, data = st, weights = sw, vcov = ~dc_cl, notes = FALSE, nthreads = 1)
        ct <- fixest::coeftable(m)
        me <- fixest::feols(f_event, data = st, weights = sw, vcov = ~dc_cl, notes = FALSE, nthreads = 1)
        ce <- fixest::coeftable(me)
        ce <- ce[grepl("^dc_erel::", rownames(ce)), , drop = FALSE]
        e <- as.numeric(sub("^dc_erel::(-?[0-9]+):dc_trt$", "\\1", rownames(ce)))
        dc_out(method, att = ct["dc_trtpost", 1], se = ct["dc_trtpost", 2],
               dynamic = dc_dyn(e, ce[, 1], ce[, 2], min_e, max_e),
               control = sprintf("Clean controls per cohort (%d sub-experiments)", length(stacks)),
               df = fixest::degrees_freedom(m, "t"), n_obs = stats::nobs(m), n_cl = length(unique(st$dc_cl)))
    }, error = function(e) dc_fail(method, e))
}
