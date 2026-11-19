# Core layer: shared runner for the single-estimator analyses (TWFE, Callaway-Sant'Anna, Sun-Abraham,
# Borusyak-Jaravel-Spiess, Gardner, stacked DiD). Each analysis class in R/did<key>.b.R is a thin
# wrapper around dc_single_run(); the estimators themselves are the dc_* functions of
# core_estimators.R, the same ones "Compare DiD Estimators" uses, so results agree exactly.

dc_single_notes <- function(key, o) {
    switch(key,
        twfe = paste0("<b>TWFE</b>: regression of the outcome on the treatment indicator with unit and period fixed effects. ",
            "With staggered adoption and heterogeneous or dynamic effects the coefficient is a weighted average that can put negative weights ",
            "on some units, so it can differ from the heterogeneity-robust estimators even in sign. The event-study version uses leads and lags ",
            "with e = -1 omitted. Use 'Bacon, Pre-trend and HonestDiD' to see the weights, or 'Compare DiD Estimators' to compare with robust estimators."),
        cs = paste0("<b>Callaway &amp; Sant'Anna</b>: group-time effects ATT(g,t) compared with ",
            if (identical(o$csControl, "notyettreated")) "not-yet-treated units" else "never-treated units",
            "; the overall ATT averages post-treatment ATT(g,t) weighted by cohort size. Universal base period (g-1). ",
            "On unbalanced panels it uses only units observed in every period. Covariates enter the propensity-score and outcome models ",
            "(method chosen above). The multiplier bootstrap is used when selected or when clustering above the unit level."),
        sa = paste0("<b>Sun &amp; Abraham</b>: interaction-weighted estimator with cohort-by-event-time effects and never-treated units as controls; ",
            "the overall ATT averages them by cohort shares. Covariates that vary within units enter linearly."),
        bjs = paste0("<b>Borusyak, Jaravel &amp; Spiess</b>: unit and period effects are fitted on untreated observations only and used to impute each ",
            "treated observation's counterfactual; the effect is the average treated-minus-imputed gap. Pre-treatment coefficients come from BJS's ",
            "pre-trend regression on untreated observations only, so they are a genuine test of the assumption behind the imputation. ",
            "Sampling weights are not available for this estimator."),
        gardner = paste0("<b>Gardner two-stage</b>: the same imputation as BJS, done as two regressions (did2s). Point estimates equal BJS; ",
            "standard errors differ because of the variance estimator. Because the first stage also uses treated units' pre-periods, its pre-treatment ",
            "coefficients are close to zero by construction and are not a pre-trend test."),
        etwfe = paste0("<b>Extended TWFE (Wooldridge)</b>: one effect for every treated cohort-by-period cell, estimated in a regression with cohort and period effects ",
            "(linear model) or with cohort and period dummies (Poisson and logit, no fixed effects); covariates enter demeaned within cohort and interacted with the cell effects, ",
            "and also with cohort and period. Effects are the average change of the predicted outcome when the treatment is switched off, on the scale of the outcome, averaged ",
            "over the treated observations (overall, by cohort, by calendar time and by event time), with delta-method standard errors. With linear models, no covariates and not-yet-treated ",
            "controls it equals the imputation estimators; with covariates it differs from them, and Poisson and logit models suit counts and 0/1 outcomes. With never-treated controls ",
            "pre-treatment effects are also estimated."),
        stacked = paste0("<b>Stacked DiD</b>: one sub-experiment per treatment cohort, using that cohort and clean controls (never treated, or not yet ",
            "treated by the end of the event window), stacked with unit-by-experiment and period-by-experiment fixed effects. The aggregate weights ",
            "cohorts by their size and by how many controls each uses, so it need not equal the other overall ATTs."))
}

dc_single_run <- function(self, key) {
    o <- self$options
    info <- self$results$info
    dc_hide_empty(self, c("dep", "id", "time", "treat"))
    self$results$guide$setContent(dc_guide_html(paste0("did", key)))
    self$results$notes$setContent(dc_method_html(paste0("did", key), o))

    if (is.null(o$dep) || is.null(o$id) || is.null(o$time) || is.null(o$treat)) {
        info$setContent(paste0(
            "<p>Select an outcome, a unit ID, a time period and a treatment ",
            "variable to estimate difference-in-differences effects.</p>"))
        return()
    }
    if (o$maxE < 0 || o$minE > -1)
        jmvcore::reject("The event window must include at least one pre-period (e = -1) and one post-period (e = 0).", code = "")

    wvar <- if (key %in% c("bjs", "etwfe")) NULL else o$weights
    d <- tryCatch(
        dc_prepare(self$data, o$dep, o$id, o$time, o$treat, o$treatType, covs = o$covs,
                   cluster = o$cluster, weights = wvar),
        error = function(e) e)
    if (inherits(d, "error"))
        jmvcore::reject(conditionMessage(d), code = "")
    info$setContent(dc_info_html(d$info))

    is_cs <- key == "cs"
    extras <- is_cs && (isTRUE(o$showGT) || isTRUE(o$aggGroup) || isTRUE(o$aggCalendar))
    r <- switch(key,
        twfe    = dc_twfe(d, o$minE, o$maxE),
        cs      = dc_cs(d, o$minE, o$maxE, o$csControl, o$csMethod,
                        boot = list(on = o$bootstrap, reps = o$bootReps, seed = o$seed), extras = extras),
        sa      = dc_sa(d, o$minE, o$maxE),
        bjs     = dc_bjs(d, o$minE, o$maxE),
        gardner = dc_gardner(d, o$minE, o$maxE),
        stacked = dc_stacked(d, o$minE, o$maxE),
        etwfe   = dc_etwfe(d, o$minE, o$maxE, o$etControl, o$etFamily))

    tab <- self$results$summary
    if (r$ok) {
        inf <- dc_inference(r$att, r$se, r$df, o$ciLevel)
        tab$addRow(rowKey = key, values = list(
            method = r$method, att = r$att, se = r$se, lower = inf$lo, upper = inf$hi, p = inf$p,
            control = r$control, nobs = r$n_obs, ncl = r$n_cl, note = ""))
    } else {
        tab$addRow(rowKey = key, values = list(method = r$method, note = r$message))
    }

    if (r$ok) {
        dtab <- self$results$dynamic
        dy <- r$dynamic
        crit <- dc_inference(0, 1, r$df, o$ciLevel)$crit
        for (i in seq_len(nrow(dy))) {
            vals <- list(e = dy$e[i], est = dy$est[i])
            if (!is.na(dy$se[i])) {
                inf <- dc_inference(dy$est[i], dy$se[i], r$df, o$ciLevel)
                vals <- c(vals, list(se = dy$se[i], lower = inf$lo, upper = inf$hi, p = inf$p))
            }
            dtab$addRow(rowKey = dy$e[i], values = vals)
        }
        long <- data.frame(method = r$method, e = dy$e, est = dy$est,
                           lo = dy$est - crit * dy$se, hi = dy$est + crit * dy$se, stringsAsFactors = FALSE)
        self$results$plot$setState(long)
    }

    if ((is_cs || key == "etwfe") && r$ok && !is.null(r$extra)) {
        z <- stats::qnorm(1 - (1 - o$ciLevel / 100) / 2)
        fill <- function(tab, df, k) {
            for (i in seq_len(nrow(df)))
                tab$addRow(rowKey = i, values = list(
                    g = as.character(dc_orig_time(d$info, df[[k]][i])), t = as.character(dc_orig_time(d$info, df[[k]][i])),
                    est = df$est[i], se = df$se[i],
                    lower = df$est[i] - z * df$se[i], upper = df$est[i] + z * df$se[i]))
        }
        if (isTRUE(o$aggGroup))    fill(self$results$groupTab, r$extra$group, "g")
        if (isTRUE(o$aggCalendar)) fill(self$results$calTab, r$extra$calendar, "t")
        if (is_cs && isTRUE(o$showGT)) {
            gt <- r$extra$gt
            gt$g <- dc_orig_time(d$info, gt$g); gt$t <- dc_orig_time(d$info, gt$t)
            self$results$gtPlot$setState(gt)
        }
    } else if (is_cs && !r$ok) {
        msg <- paste("Callaway & Sant'Anna did not run:", r$message)
        if (isTRUE(o$aggGroup))    self$results$groupTab$setNote("err", msg)
        if (isTRUE(o$aggCalendar)) self$results$calTab$setNote("err", msg)
    }

    if (isTRUE(o$showCode) && key == "etwfe") {
        self$results$code$setContent(dc_code_etwfe(list(dep = o$dep, id = o$id, time = o$time, treat = o$treat, treatType = o$treatType, covs = o$covs, cluster = o$cluster,
                                                        control = o$etControl, family = o$etFamily, minE = o$minE, maxE = o$maxE), d$info))
    } else if (isTRUE(o$showCode)) {
        flags <- stats::setNames(as.list(names(c(twfe = 1, cs = 1, sa = 1, bjs = 1, gardner = 1, stacked = 1)) == key),
                                 c("twfe", "cs", "sa", "bjs", "gardner", "stacked"))
        opts <- c(list(dep = o$dep, id = o$id, time = o$time, treat = o$treat, treatType = o$treatType,
                       covs = o$covs, cluster = o$cluster, weights = wvar, minE = o$minE, maxE = o$maxE,
                       csControl = if (is_cs) o$csControl else "nevertreated",
                       csMethod = if (is_cs) o$csMethod else "dr",
                       boot = is_cs && isTRUE(o$bootstrap),
                       bootReps = if (is_cs) o$bootReps else 999L,
                       seed = if (is_cs) o$seed else 1234L), flags)
        self$results$code$setContent(dc_code_compare(opts, d$info))
    }

    run <- paste0("<p>", r$method, ": ", if (r$ok) "estimated" else paste("did not run (", r$message, ")"),
        " on ", d$info$n_units, " units, standard errors clustered at the ",
        if (length(o$cluster)) "level of the cluster variable" else "unit level", ".</p>")
    self$results$notes$setContent(dc_run_notes(paste0("did", key), o, run))
}

dc_single_plot <- function(image, ggtheme, o) {
    state <- image$state
    if (is.null(state))
        return(FALSE)
    print(dc_plot_event(state, ggtheme, o$plotPalette, o$plotCI, 0.5))
    TRUE
}
