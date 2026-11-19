
didcompareClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didcompareClass",
    inherit = didcompareBase,
    private = list(

        .run = function() {
            self$results$guide$setContent(dc_guide_html("didcompare"))
            dc_hide_empty(self, c("dep", "id", "time", "treat"))
            self$results$notes$setContent(dc_method_html("didcompare", self$options))
            o <- self$options
            info <- self$results$info

            if (is.null(o$dep) || is.null(o$id) || is.null(o$time) || is.null(o$treat)) {
                info$setContent(paste0(
                    "<p>Select an outcome, a unit ID, a time period and a treatment ",
                    "variable to compare difference-in-differences estimators.</p>"))
                return()
            }
            if (!(o$twfe || o$cs || o$sa || o$bjs || o$gardner || o$stacked)) {
                info$setContent("<p>Select at least one estimator.</p>")
                return()
            }
            if (o$maxE < 0 || o$minE > -1)
                jmvcore::reject("The event window must include at least one pre-period (e = -1) and one post-period (e = 0).", code = "")

            d <- tryCatch(
                dc_prepare(self$data, o$dep, o$id, o$time, o$treat, o$treatType, covs = o$covs,
                           cluster = o$cluster, weights = o$weights),
                error = function(e) e)
            if (inherits(d, "error"))
                jmvcore::reject(conditionMessage(d), code = "")

            info$setContent(dc_info_html(d$info))

            res <- list()
            if (o$twfe) res$twfe <- dc_twfe(d, o$minE, o$maxE)
            boot <- list(on = o$bootstrap, reps = o$bootReps, seed = o$seed)
            if (o$cs)   res$cs   <- dc_cs(d, o$minE, o$maxE, o$csControl, o$csMethod, boot = boot,
                                          extras = o$showGT || o$aggGroup || o$aggCalendar)
            if (o$sa)   res$sa   <- dc_sa(d, o$minE, o$maxE)
            if (o$bjs)  res$bjs  <- dc_bjs(d, o$minE, o$maxE)
            if (o$gardner) res$gardner <- dc_gardner(d, o$minE, o$maxE)
            if (o$stacked) res$stacked <- dc_stacked(d, o$minE, o$maxE)

            z <- stats::qnorm(1 - (1 - o$ciLevel / 100) / 2)

            # ---- overall ATT table ----
            tab <- self$results$summary
            for (key in names(res)) {
                r <- res[[key]]
                if (r$ok) {
                    inf <- dc_inference(r$att, r$se, r$df, o$ciLevel)
                    tab$addRow(rowKey = key, values = list(
                        method = r$method, att = r$att, se = r$se,
                        lower = inf$lo, upper = inf$hi, p = inf$p,
                        control = r$control, nobs = r$n_obs, ncl = r$n_cl, note = ""))
                } else {
                    tab$addRow(rowKey = key, values = list(
                        method = r$method, note = r$message))
                }
            }

            # ---- dynamic table ----
            dtab <- self$results$dynamic
            for (e in seq(o$minE, o$maxE)) {
                vals <- list(e = e)
                for (key in names(res)) {
                    r <- res[[key]]
                    if (!r$ok) next
                    i <- match(e, r$dynamic$e)
                    if (is.na(i)) next
                    vals[[paste0(key, "_est")]] <- r$dynamic$est[i]
                    if (!is.na(r$dynamic$se[i]))
                        vals[[paste0(key, "_se")]] <- r$dynamic$se[i]
                }
                dtab$addRow(rowKey = e, values = vals)
            }

            # ---- plot data ----
            long <- do.call(rbind, lapply(res, function(r) {
                if (!r$ok) return(NULL)
                crit <- dc_inference(0, 1, r$df, o$ciLevel)$crit     # each estimator's own reference distribution
                data.frame(method = r$method, e = r$dynamic$e, est = r$dynamic$est,
                           lo = r$dynamic$est - crit * r$dynamic$se,
                           hi = r$dynamic$est + crit * r$dynamic$se,
                           stringsAsFactors = FALSE)
            }))
            if (!is.null(long))
                self$results$plot$setState(long)

            # ---- forest plot ----
            forest <- do.call(rbind, lapply(res, function(r) {
                if (!r$ok) return(NULL)
                inf <- dc_inference(r$att, r$se, r$df, o$ciLevel)
                data.frame(method = r$method, est = r$att, lo = inf$lo, hi = inf$hi, stringsAsFactors = FALSE)
            }))
            if (!is.null(forest))
                self$results$forest$setState(forest)

            # ---- Callaway & Sant'Anna aggregations ----
            cs <- res$cs
            if (!is.null(cs) && cs$ok && !is.null(cs$extra)) {
                fill <- function(tab, df, key) {
                    for (i in seq_len(nrow(df)))
                        tab$addRow(rowKey = i, values = list(
                            g = as.character(dc_orig_time(d$info, df[[key]][i])), t = as.character(dc_orig_time(d$info, df[[key]][i])),
                            est = df$est[i], se = df$se[i],
                            lower = df$est[i] - z * df$se[i], upper = df$est[i] + z * df$se[i]))
                }
                if (o$aggGroup)    fill(self$results$groupTab, cs$extra$group, "g")
                if (o$aggCalendar) fill(self$results$calTab, cs$extra$calendar, "t")
                if (o$showGT) {
                    gt <- cs$extra$gt
                    gt$g <- dc_orig_time(d$info, gt$g); gt$t <- dc_orig_time(d$info, gt$t)
                    self$results$gtPlot$setState(gt)
                }
            } else if (!is.null(cs) && !cs$ok) {
                msg <- paste("Callaway & Sant'Anna did not run:", cs$message)
                if (o$aggGroup)    self$results$groupTab$setNote("err", msg)
                if (o$aggCalendar) self$results$calTab$setNote("err", msg)
            }

            # ---- reproducible code ----
            if (o$showCode) {
                opts <- list(dep = o$dep, id = o$id, time = o$time, treat = o$treat, treatType = o$treatType,
                             covs = o$covs, cluster = o$cluster, weights = o$weights, minE = o$minE, maxE = o$maxE,
                             twfe = o$twfe, cs = o$cs, sa = o$sa, bjs = o$bjs, gardner = o$gardner,
                             stacked = o$stacked, csControl = o$csControl, csMethod = o$csMethod,
                             boot = o$bootstrap, bootReps = o$bootReps, seed = o$seed)
                self$results$code$setContent(dc_code_compare(opts, d$info))
            }

            self$results$notes$setContent(dc_run_notes("didcompare", self$options, private$.notesHtml(o, res)))
        },

        .plotForest = function(image, ggtheme, theme, ...) {
            state <- image$state
            if (is.null(state))
                return(FALSE)
            print(dc_plot_forest(state, ggtheme, self$options$plotPalette))
            TRUE
        },

        .plotGT = function(image, ggtheme, theme, ...) {
            state <- image$state
            if (is.null(state))
                return(FALSE)
            print(dc_plot_gt(state, ggtheme, self$options$plotPalette))
            TRUE
        },

        .plot = function(image, ggtheme, theme, ...) {
            state <- image$state
            if (is.null(state))
                return(FALSE)
            print(dc_plot_event(state, ggtheme, self$options$plotPalette, self$options$plotCI, self$options$plotDodge))
            TRUE
        },

        .notesHtml = function(o, res) {
            n <- character()
            if (!is.null(res$twfe))
                n <- c(n, paste0("<b>TWFE</b>: regression of the outcome on the treatment indicator with unit and ",
                    "period fixed effects. With staggered adoption and heterogeneous or dynamic effects the ",
                    "coefficient is a weighted average that can put negative weights on some units, so it can differ ",
                    "from the other estimators even in sign. The event-study version uses leads and lags with e = -1 omitted."))
            if (!is.null(res$cs))
                n <- c(n, paste0("<b>Callaway &amp; Sant'Anna</b>: group-time effects ATT(g,t) compared with ",
                    if (o$csControl == "nevertreated") "never-treated units" else "not-yet-treated units",
                    "; the overall ATT averages post-treatment ATT(g,t) weighted by cohort size. ",
                    "Universal base period (g-1). On unbalanced panels it uses only units observed in every period."))
            if (!is.null(res$sa))
                n <- c(n, paste0("<b>Sun &amp; Abraham</b>: interaction-weighted estimator with cohort-by-event-time ",
                    "effects and never-treated units as controls; the overall ATT averages them by cohort shares."))
            if (!is.null(res$bjs))
                n <- c(n, paste0("<b>Borusyak, Jaravel &amp; Spiess</b>: unit and period effects are fitted on untreated ",
                    "observations only and used to impute each treated observation's counterfactual; the effect is the average ",
                    "treated-minus-imputed gap. The overall ATT averages over all treated observations. Pre-treatment coefficients come from ",
                    "BJS's pre-trend regression on untreated observations only, so they are a genuine test of the assumption behind the imputation."))
            if (!is.null(res$gardner))
                n <- c(n, paste0("<b>Gardner two-stage</b>: the same imputation as BJS, done as two regressions (did2s). ",
                    "Point estimates equal BJS; standard errors differ because of the variance estimator. Because the first stage also ",
                    "uses treated units' pre-periods, its pre-treatment coefficients are close to zero by construction and are not a pre-trend test."))
            if (!is.null(res$stacked))
                n <- c(n, paste0("<b>Stacked DiD</b>: one sub-experiment per treatment cohort, using that cohort and clean controls ",
                    "(never treated, or not yet treated by the end of the event window), stacked with unit-by-experiment and ",
                    "period-by-experiment fixed effects. The aggregate weights cohorts by their size and by how many controls each uses, ",
                    "so it need not equal the other overall ATTs."))
            if (!is.null(res$bjs) || !is.null(res$gardner))
                n <- c(n, paste0("Wooldridge's extended TWFE is not listed separately: with never-treated controls its post-treatment ",
                    "effects equal the imputation estimates above, and with lead effects added it equals Sun &amp; Abraham."))
            if (length(o$cluster) || o$bootstrap)
                n <- c(n, paste0("<b>Standard errors</b>: clustered at the ",
                    if (length(o$cluster)) "level of the cluster variable" else "unit level",
                    ". Callaway &amp; Sant'Anna uses its multiplier bootstrap (", o$bootReps, " draws, seed ", o$seed,
                    ") when bootstrap is selected or when clustering above the unit level, because the did package ",
                    "only supports such clustering through the bootstrap; the other estimators cluster analytically."))
            if (length(o$covs))
                n <- c(n, paste0("<b>Covariates</b>: Callaway &amp; Sant'Anna uses them in the propensity-score and outcome models ",
                    "(method chosen above) and so can handle trends that depend on baseline characteristics. The other estimators enter ",
                    "covariates linearly in the regression (or first stage) next to the fixed effects, so only covariates that change over ",
                    "time within a unit matter for them."))
            n <- c(n, paste0("<b>Inference</b>: TWFE, Sun-Abraham and stacked DiD use t(G-1) critical values and p-values with G clusters, as fixest and Stata's ",
                "reghdfe do; Callaway-Sant'Anna, Borusyak-Jaravel-Spiess and Gardner use their packages' normal approximation. With few clusters ",
                "(see the warning above) all of these can be too optimistic."))
            n <- c(n, paste0("Standard errors are clustered at the unit level. Differences between estimators can ",
                "reflect different weighting of the same effects, different control groups, or bias from the ",
                "TWFE specification, not only sampling error."))
            paste0("<p>", paste(n, collapse = "</p><p>"), "</p>")
        })
)
