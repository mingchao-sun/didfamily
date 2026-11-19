# Core layer: classic two-group, two-period DiD. No jamovi dependencies.
#
# `group`: 0/1 variable (1 = treated group); `period`: variable with exactly two values, the
# later / larger / last level being the post period; `id`: optional unit id (panel data with
# both periods for every unit; without it the data are treated as repeated cross-sections).

dc_prepare22 <- function(data, y, group, period, id = NULL, covs = NULL) {
    pv <- data[[period]]
    lev <- if (is.factor(pv)) levels(droplevels(pv)) else sort(unique(pv[!is.na(pv)]))
    if (length(lev) != 2)
        stop(sprintf("The period variable must have exactly two values (before and after); it has %d.", length(lev)))
    df <- data.frame(
        dc_y = dc_num(data[[y]]),
        dc_g = dc_num(data[[group]]),
        dc_post = as.numeric(as.character(pv) == as.character(lev[2])),
        stringsAsFactors = FALSE)
    df$dc_post[is.na(pv)] <- NA
    b <- dc_first_bad(data[[y]], df$dc_y)
    if (!is.null(b))
        stop(sprintf("The outcome variable must be numeric: %s values such as '%s' are not numbers.", b[1], b[2]))
    b <- dc_first_bad(data[[group]], df$dc_g)
    if (!is.null(b))
        stop(sprintf("The treated-group variable must be numeric 0/1: %s values such as '%s' are not numbers. Recode it to 0 and 1.", b[1], b[2]))
    if (!is.null(id))
        df$dc_id <- as.character(data[[id]])
    cb <- dc_build_covs(data, covs)
    for (cd in names(cb$cols)) df[[cd]] <- cb$cols[[cd]]

    n_input <- nrow(df)
    df <- df[stats::complete.cases(df), , drop = FALSE]
    n_missing <- n_input - nrow(df)
    if (!all(df$dc_g %in% c(0, 1)))
        stop("The group variable must contain only 0 (comparison group) and 1 (treated group).")
    panel <- !is.null(id)
    if (panel) {
        df$dc_id <- match(df$dc_id, unique(df$dc_id))
        if (anyDuplicated(df[, c("dc_id", "dc_post")]))
            stop("Some units have more than one row in the same period.")
        cnt <- table(df$dc_id)
        if (any(cnt != 2))
            stop("With a unit ID every unit must be observed in both periods (balanced two-period panel). ",
                 "Remove the unit ID to treat the data as repeated cross-sections.")
        if (any(stats::ave(df$dc_g, df$dc_id, FUN = function(v) length(unique(v))) > 1))
            stop("The group variable must be constant within each unit.")
        df <- df[order(df$dc_id, df$dc_post), , drop = FALSE]
    }
    cells <- table(df$dc_g, df$dc_post)
    if (any(dim(cells) < 2) || any(cells == 0))
        stop("Each combination of group and period needs at least one observation.")
    covmap <- cb$map
    list(data = df, panel = panel, labels = as.character(lev), covs = covmap,
         info = list(n_obs = nrow(df), n_units = if (panel) length(unique(df$dc_id)) else NA_integer_,
                     n_missing = n_missing, cells = cells))
}

# cell means and sizes: one row per group x period
dc_means22 <- function(d) {
    df <- d$data
    out <- expand.grid(group = c(0, 1), post = c(0, 1))
    out$mean <- mapply(function(g, p) mean(df$dc_y[df$dc_g == g & df$dc_post == p]), out$group, out$post)
    out$n <- mapply(function(g, p) sum(df$dc_g == g & df$dc_post == p), out$group, out$post)
    out
}

# `methods`: any of "ols", "dr", "drtrad", "ipw", "or". Each result: list(method, ok, message, att, se)
dc_did22 <- function(d, methods = c("ols", "dr", "drtrad", "ipw", "or")) {
    df <- d$data
    codes <- if (nrow(d$covs)) d$covs$code else character()
    cx <- if (length(codes)) str2lang(paste(codes, collapse = " + ")) else NULL
    xf <- if (is.null(cx)) ~1 else eval(bquote(~ .(cx)))
    if (!d$panel) df$dc_id <- seq_len(nrow(df))     # DRDID wants an id column even for cross-sections
    fail <- function(m, e) list(method = m, ok = FALSE, message = conditionMessage(e), att = NA_real_, se = NA_real_)
    labels <- c(ols = "Regression DiD (OLS)", dr = "Doubly robust (improved)", drtrad = "Doubly robust (traditional)",
                ipw = "Inverse probability weighting", or = "Outcome regression")
    out <- list()
    for (m in methods) {
        out[[m]] <- tryCatch({
            if (m == "ols") {
                f <- if (is.null(cx)) dc_y ~ dc_g * dc_post else eval(bquote(dc_y ~ dc_g * dc_post + .(cx)))
                fit <- fixest::feols(f, data = df, vcov = if (d$panel) ~dc_id else "hetero",
                                     notes = FALSE, nthreads = 1)
                ct <- fixest::coeftable(fit)
                list(method = labels[[m]], ok = TRUE, message = "", att = ct["dc_g:dc_post", 1], se = ct["dc_g:dc_post", 2])
            } else {
                fn <- switch(m, dr = DRDID::drdid, drtrad = DRDID::drdid, ipw = DRDID::ipwdid, or = DRDID::ordid)
                args <- list(yname = "dc_y", tname = "dc_post", idname = "dc_id", dname = "dc_g",
                             xformla = xf, data = df, panel = d$panel)
                if (m == "dr") args$estMethod <- "imp"
                if (m == "drtrad") args$estMethod <- "trad"
                fit <- do.call(fn, args)
                list(method = labels[[m]], ok = TRUE, message = "", att = fit$ATT, se = fit$se)
            }
        }, error = function(e) fail(labels[[m]], e))
    }
    out
}
