
did22Class <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "did22Class",
    inherit = did22Base,
    private = list(

        .run = function() {
            self$results$guide$setContent(dc_guide_html("did22"))
            dc_hide_empty(self, c("dep", "group", "period"))
            self$results$notes$setContent(dc_method_html("did22", self$options))
            o <- self$options
            info <- self$results$info

            if (is.null(o$dep) || is.null(o$group) || is.null(o$period)) {
                info$setContent(paste0(
                    "<p>Select an outcome, a 0/1 variable marking the treated group, and a variable with the two periods ",
                    "(before and after). Add a unit ID if the same units are observed in both periods.</p>"))
                return()
            }
            methods <- c(if (o$ols) "ols", if (o$dr) "dr", if (o$drtrad) "drtrad", if (o$ipw) "ipw", if (o$or) "or")
            if (!length(methods)) {
                info$setContent("<p>Select at least one estimator.</p>")
                return()
            }

            d <- tryCatch(
                dc_prepare22(self$data, o$dep, o$group, o$period, id = o$id, covs = o$covs),
                error = function(e) e)
            if (inherits(d, "error"))
                jmvcore::reject(conditionMessage(d), code = "")

            msg <- c(sprintf("<b>%d</b> observations (%s). Before = <b>%s</b>, after = <b>%s</b>.",
                             d$info$n_obs,
                             if (d$panel) sprintf("balanced panel of %d units", d$info$n_units) else "repeated cross-sections",
                             d$labels[1], d$labels[2]))
            if (nrow(d$covs))
                msg <- c(msg, sprintf("Covariates: %s.", paste(d$covs$label, collapse = ", ")))
            if (d$info$n_missing > 0)
                msg <- c(msg, sprintf("Dropped %d rows with missing values.", d$info$n_missing))
            info$setContent(paste0("<p>", paste(msg, collapse = "<br>"), "</p>"))

            # ---- means ----
            m <- dc_means22(d)
            tab <- self$results$means
            for (g in c(0, 1)) {
                pre <- m$mean[m$group == g & m$post == 0]; post <- m$mean[m$group == g & m$post == 1]
                tab$addRow(rowKey = g, values = list(
                    group = if (g == 1) "Treated" else "Comparison", pre = pre, post = post,
                    change = post - pre, n = sum(m$n[m$group == g])))
            }
            att_raw <- (m$mean[m$group == 1 & m$post == 1] - m$mean[m$group == 1 & m$post == 0]) -
                       (m$mean[m$group == 0 & m$post == 1] - m$mean[m$group == 0 & m$post == 0])
            tab$addRow(rowKey = "did", values = list(group = "Difference in changes", change = att_raw))

            # ---- estimates ----
            res <- dc_did22(d, methods)
            z <- stats::qnorm(1 - (1 - o$ciLevel / 100) / 2)
            stab <- self$results$summary
            for (key in names(res)) {
                r <- res[[key]]
                if (r$ok) {
                    stab$addRow(rowKey = key, values = list(
                        method = r$method, att = r$att, se = r$se,
                        lower = r$att - z * r$se, upper = r$att + z * r$se,
                        p = 2 * stats::pnorm(-abs(r$att / r$se)), note = ""))
                } else {
                    stab$addRow(rowKey = key, values = list(method = r$method, note = r$message))
                }
            }

            self$results$plot$setState(list(means = m, labels = d$labels))
            self$results$notes$setContent(dc_run_notes("did22", self$options, private$.notesHtml(o, d)))
        },

        .plot = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            print(dc_plot_22(s$means, s$labels, ggtheme, self$options$plotPalette))
            TRUE
        },

        .notesHtml = function(o, d) {
            n <- c(
                paste0("<b>Regression DiD</b>: OLS of the outcome on group, period and their interaction (plus covariates entered ",
                    "linearly); the interaction is the ATT and equals the difference in changes of the group means when there are no covariates."),
                paste0("<b>Doubly robust, IPW and outcome regression</b> (Sant'Anna &amp; Zhao 2020) use covariates to make the groups comparable ",
                    "before comparing their changes: IPW reweights the comparison group by a propensity score, outcome regression predicts the ",
                    "treated group's change from the comparison group, and the doubly robust estimators combine both and stay consistent if either model is right. ",
                    "Without covariates all estimators coincide."),
                paste0("The identifying assumption is that, given the covariates, the treated group would have changed like the comparison group (parallel trends). ",
                    "With only two periods it cannot be tested from this data; the dashed line in the plot shows the counterfactual it implies."),
                if (d$panel) "Standard errors are clustered by unit." else "Standard errors are heteroskedasticity-robust (repeated cross-sections).")
            paste0("<p>", paste(n, collapse = "</p><p>"), "</p>")
        })
)
