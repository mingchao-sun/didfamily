
dddcompareClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "dddcompareClass",
    inherit = dddcompareBase,
    private = list(

        .run = function() {
            self$results$guide$setContent(dc_guide_html("dddcompare"))
            dc_hide_empty(self, c("dep", "id", "time", "treat", "part"))
            self$results$notes$setContent(dc_method_html("dddcompare", self$options))
            o <- self$options
            info <- self$results$info

            if (is.null(o$dep) || is.null(o$id) || is.null(o$time) || is.null(o$treat) || is.null(o$part)) {
                info$setContent(paste0(
                    "<p>Select an outcome, a unit ID, a time period, the first treated period of the unit's group ",
                    "(0 for never treated) and a 0/1 eligibility variable to estimate triple differences.</p>"))
                return()
            }
            if (o$maxE < 0 || o$minE > -1)
                jmvcore::reject("The event window must include at least one pre-period (e = -1) and one post-period (e = 0).", code = "")

            d <- tryCatch(
                dc_prepare(self$data, o$dep, o$id, o$time, o$treat, "cohort",
                           covs = o$covs, part = o$part),
                error = function(e) e)
            if (inherits(d, "error"))
                jmvcore::reject(conditionMessage(d), code = "")
            n_el <- length(unique(d$data$dc_id[d$data$dc_p == 1]))
            info$setContent(paste0(dc_info_html(d$info),
                sprintf("<p>Eligible units: <b>%d</b>; ineligible units: <b>%d</b>.</p>",
                        n_el, d$info$n_units - n_el)))
            if (n_el == 0 || n_el == d$info$n_units)
                jmvcore::reject("The eligibility variable must contain both eligible (1) and ineligible (0) units.", code = "")

            res <- list(
                ddd  = dc_ddd(d, o$minE, o$maxE, o$csControl, o$csMethod),
                elig = dc_did_part(d, "eligible", o$minE, o$maxE, o$csControl, o$csMethod),
                inel = dc_did_part(d, "ineligible", o$minE, o$maxE, o$csControl, o$csMethod))
            private$.checkpoint()

            z <- stats::qnorm(1 - (1 - o$ciLevel / 100) / 2)
            tab <- self$results$summary
            for (key in names(res)) {
                r <- res[[key]]
                if (r$ok) {
                    tab$addRow(rowKey = key, values = list(
                        method = r$method, att = r$att, se = r$se,
                        lower = r$att - z * r$se, upper = r$att + z * r$se,
                        p = 2 * stats::pnorm(-abs(r$att / r$se)),
                        control = r$control, note = ""))
                } else {
                    tab$addRow(rowKey = key, values = list(method = r$method, note = r$message))
                }
            }

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

            long <- do.call(rbind, lapply(res, function(r) {
                if (!r$ok) return(NULL)
                data.frame(method = r$method, e = r$dynamic$e, est = r$dynamic$est,
                           lo = r$dynamic$est - z * r$dynamic$se,
                           hi = r$dynamic$est + z * r$dynamic$se,
                           stringsAsFactors = FALSE)
            }))
            if (!is.null(long))
                self$results$plot$setState(long)

            self$results$notes$setContent(dc_run_notes("dddcompare", self$options, private$.notesHtml(o)))
        },

        .plot = function(image, ggtheme, theme, ...) {
            state <- image$state
            if (is.null(state))
                return(FALSE)
            print(dc_plot_event(state, ggtheme, self$options$plotPalette, self$options$plotCI, self$options$plotDodge))
            TRUE
        },

        .notesHtml = function(o) {
            n <- c(
                paste0("<b>Triple difference (DDD)</b>: compares how the outcome of <i>eligible</i> units changes after the policy ",
                    "with how the outcome of <i>ineligible</i> units changes, in treated relative to control groups (Ortiz-Villavicencio ",
                    "&amp; Sant'Anna). It stays valid if a shock hits the whole treated group, as long as that shock would have affected ",
                    "eligible and ineligible units equally."),
                paste0("<b>DiD, eligible only</b> ignores the ineligible units. If treated groups also experience a group-wide shock, this ",
                    "estimate absorbs it. <b>DiD, ineligible only</b> should be zero if the policy does not affect them; a non-zero ",
                    "estimate points to exactly such a shock, and is the amount by which the eligible-only DiD is off."),
                paste0("Group-time effects use the ", if (o$csControl == "nevertreated") "never-treated" else "not-yet-treated",
                    " groups as controls with a universal base period (g-1). Covariates enter the propensity-score and outcome models ",
                    "(method chosen above). Standard errors are analytical (influence-function based)."))
            paste0("<p>", paste(n, collapse = "</p><p>"), "</p>")
        })
)
