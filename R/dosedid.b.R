
dosedidClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "dosedidClass",
    inherit = dosedidBase,
    private = list(

        .run = function() {
            self$results$guide$setContent(dc_guide_html("dosedid"))
            dc_hide_empty(self, c("dep", "id", "time", "treat", "dose"))
            self$results$notes$setContent(dc_method_html("dosedid", self$options))
            o <- self$options
            info <- self$results$info

            if (is.null(o$dep) || is.null(o$id) || is.null(o$time) || is.null(o$treat) || is.null(o$dose)) {
                info$setContent(paste0(
                    "<p>Select an outcome, a unit ID, a time period, the first treated period (0 for never treated) ",
                    "and a dose variable (the amount of treatment a unit receives once treated, 0 if untreated).</p>"))
                return()
            }
            d <- tryCatch(
                dc_prepare(self$data, o$dep, o$id, o$time, o$treat, "cohort", dose = o$dose),
                error = function(e) e)
            if (inherits(d, "error"))
                jmvcore::reject(conditionMessage(d), code = "")
            tdose <- d$data$dc_dose[!duplicated(d$data$dc_id) & d$data$dc_g > 0]
            info$setContent(paste0(dc_info_html(d$info), if (length(tdose))
                sprintf("<p>Dose among treated units: min %.3g, median %.3g, max %.3g.</p>",
                        min(tdose), stats::median(tdose), max(tdose))))

            r <- dc_dose(d, o$csControl, o$degree, o$knots)
            private$.checkpoint()
            tab <- self$results$summary
            if (!r$ok) {
                tab$addRow(rowKey = "err", values = list(param = "Continuous-treatment DiD", note = r$message))
                return()
            }
            z <- stats::qnorm(1 - (1 - o$ciLevel / 100) / 2)
            add <- function(key, name, est, se)
                tab$addRow(rowKey = key, values = list(
                    param = name, est = est, se = se, lower = est - z * se, upper = est + z * se,
                    p = 2 * stats::pnorm(-abs(est / se)), note = ""))
            add("att", "Average effect of the treatment (ATT)", r$att, r$att_se)
            add("acrt", "Average marginal effect of dose (ACRT)", r$acrt, r$acrt_se)

            cv <- r$curve
            show <- cv[unique(round(seq(1, nrow(cv), length.out = min(10, nrow(cv))))), ]
            ct <- self$results$curve
            for (i in seq_len(nrow(show)))
                ct$addRow(rowKey = i, values = as.list(show[i, ]))
            self$results$plot$setState(list(curve = cv, z = z))

            self$results$notes$setContent(dc_run_notes("dosedid", self$options, paste0(
                "<p><b>ATT(d)</b> is the average effect on units that received dose d, compared with what would have happened with no treatment. ",
                "<b>ACRT(d)</b> is the average marginal effect: how much the outcome changes per extra unit of dose around d. ",
                "The overall ATT and ACRT average these over the treated units' doses.</p>",
                "<p>Estimates use a B-spline dose-response of the chosen degree and number of interior knots (Callaway, Goodman-Bacon &amp; Sant'Anna 2024). ",
                "The strong version of parallel trends needed to compare effects across doses is stronger than the usual assumption for a binary treatment: ",
                "units with different doses must have shared the same untreated trend, so ATT(d) differences can also reflect selection into dose. ",
                "Standard errors are analytical and pointwise (not simultaneous bands); covariates are not used here.</p>",
                "<p>Treated units' doses start well above zero while control units sit at dose 0, so a low-degree curve can be too smooth near the smallest doses. ",
                "If the curve looks off there, raise the spline degree or add interior knots and check that the shape is stable.</p>",
                if (!is.null(r$dropped))
                    sprintf("<p>%d units with missing periods were set aside: the estimator needs units observed in every period.</p>", r$dropped))))
        },

        .plot = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            print(dc_plot_dose(s$curve, s$z, ggtheme, self$options$plotPalette))
            TRUE
        })
)
