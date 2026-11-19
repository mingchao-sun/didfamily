
didiagClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didiagClass",
    inherit = didiagBase,
    private = list(

        .run = function() {
            self$results$guide$setContent(dc_guide_html("didiag"))
            dc_hide_empty(self, c("dep", "id", "time", "treat"))
            self$results$notes$setContent(dc_method_html("didiag", self$options))
            o <- self$options
            info <- self$results$info

            if (is.null(o$dep) || is.null(o$id) || is.null(o$time) || is.null(o$treat)) {
                info$setContent(paste0(
                    "<p>Select an outcome, a unit ID, a time period and a treatment ",
                    "variable to run difference-in-differences diagnostics.</p>"))
                return()
            }
            if (!(o$bacon || o$weights || o$pretest || o$honest || o$ptpower)) {
                info$setContent("<p>Select at least one diagnostic.</p>")
                return()
            }
            if (o$maxE < 0 || o$minE > -1)
                jmvcore::reject("The event window must include at least one pre-period (e = -1) and one post-period (e = 0).", code = "")

            d <- tryCatch(
                dc_prepare(self$data, o$dep, o$id, o$time, o$treat, o$treatType),
                error = function(e) e)
            if (inherits(d, "error"))
                jmvcore::reject(conditionMessage(d), code = "")
            info$setContent(dc_info_html(d$info))

            facts <- list()
            if (o$bacon)   facts$bacon   <- private$.runBacon(d)
            if (o$weights) facts$weights <- private$.runWeights(d)
            if (o$pretest) private$.runPretest(d)
            if (o$honest)  private$.runHonest(d)
            if (o$ptpower) private$.runPtPower(d)
            self$results$notes$setContent(dc_run_notes("didiag", self$options, private$.notesHtml(o, facts)))
        },

        .runBacon = function(d) {
            tab <- self$results$bacon
            b <- dc_bacon(d)
            if (!b$ok) {
                tab$addRow(rowKey = "err", values = list(type = b$message))
                return(NULL)
            }
            for (i in seq_len(nrow(b$summary)))
                tab$addRow(rowKey = i, values = list(
                    type = b$summary$type[i], n = b$summary$n[i],
                    weight = b$summary$weight[i], est = b$summary$est[i]))
            tab$addRow(rowKey = "total", values = list(
                type = "Static TWFE estimate", n = nrow(b$detail), weight = 1, est = b$twfe))
            self$results$baconPlot$setState(list(detail = b$detail, twfe = b$twfe))
            bad <- b$summary$weight[b$summary$type == "Later vs Earlier Treated"]
            list(bad = if (length(bad)) bad else 0)
        },

        .runWeights = function(d) {
            tab <- self$results$weights
            w <- dc_twfe_weights(d)
            if (!w$ok) {
                tab$addRow(rowKey = "err", values = list(group = w$message))
                return(NULL)
            }
            tab$addRow(rowKey = "all", values = list(
                group = "All treated observations", n = w$n_treated, weight = 1,
                neg_n = w$neg_n, neg_weight = w$neg_w))
            for (i in seq_len(nrow(w$by_cohort)))
                tab$addRow(rowKey = i, values = list(
                    group = as.character(dc_orig_time(d$info, w$by_cohort$cohort[i])), n = w$by_cohort$n[i],
                    weight = w$by_cohort$weight[i], neg_n = w$by_cohort$neg_n[i],
                    neg_weight = w$by_cohort$neg_weight[i]))
            list(neg_w = w$neg_w, neg_n = w$neg_n, n = w$n_treated)
        },

        .runPretest = function(d) {
            tab <- self$results$pretest
            pt <- dc_pretests(d, self$options$minE, self$options$csControl)
            for (i in seq_len(nrow(pt))) {
                vals <- list(method = pt$method[i], note = pt$note[i])
                if (!is.na(pt$p[i])) {
                    vals$stat <- pt$stat[i]; vals$df <- pt$df[i]; vals$p <- pt$p[i]
                }
                tab$addRow(rowKey = i, values = vals)
            }
        },

        .runHonest = function(d) {
            o <- self$options
            tab <- self$results$honest
            prep <- tryCatch(
                dc_honest_prepare(d, o$honestBase, o$minE, o$maxE, o$csControl),
                error = function(e) e)
            if (inherits(prep, "error")) {
                tab$addRow(rowKey = "err", values = list(spec = "", note = conditionMessage(prep)))
                return()
            }
            alpha <- 1 - o$ciLevel / 100
            grid <- c(NA, seq(0, o$honestMbar, length.out = 5))
            labels <- c("Original (no violation allowed)",
                        sprintf("Mbar = %.2f", grid[-1]))
            labels[2] <- paste0(labels[2], " (parallel trends)")
            rows <- list()
            for (k in seq_along(grid)) {
                r <- tryCatch(dc_honest_one(prep, grid[k], o$honestTarget, alpha),
                              error = function(e) e)
                if (inherits(r, "error")) {
                    tab$addRow(rowKey = k, values = list(spec = labels[k], note = conditionMessage(r)))
                } else {
                    tab$addRow(rowKey = k, values = list(spec = labels[k], lower = r[["lb"]], upper = r[["ub"]]))
                    rows[[length(rows) + 1]] <- data.frame(
                        label = if (is.na(grid[k])) "Original" else sprintf("Mbar = %.2f", grid[k]),
                        lb = r[["lb"]], ub = r[["ub"]], original = is.na(grid[k]))
                }
                private$.checkpoint()
            }
            if (length(rows))
                self$results$honestPlot$setState(do.call(rbind, rows))
        },

        .runPtPower = function(d) {
            o <- self$options
            tab <- self$results$ptTable
            prep <- tryCatch(dc_honest_prepare(d, o$ptBase, o$minE, o$maxE, o$csControl), error = function(e) e)
            if (inherits(prep, "error")) {
                tab$addRow(rowKey = "err", values = list(item = "", note = conditionMessage(prep)))
                return()
            }
            pr <- tryCatch(dc_pretrend(prep, if (o$ptSlope != 0) o$ptSlope else NA_real_, o$ciLevel / 100), error = function(e) e)
            if (inherits(pr, "error")) {
                tab$addRow(rowKey = "err", values = list(item = "", note = conditionMessage(pr)))
                return()
            }
            k <- length(pr$e_pre)
            add <- function(i, item, value, note = "") tab$addRow(rowKey = i, values = list(item = item, value = value, note = note))
            add(1, "Largest |t| among the pre-treatment effects", pr$max_t, sprintf("critical value %s: the test %s", dc_fmt3(pr$crit), if (pr$passes) "passes (no pre-treatment effect is significant)" else "rejects (at least one pre-treatment effect is significant)"))
            add(2, "Joint Wald test of the pre-treatment effects, p", pr$wald_p, sprintf("chi-square %s with %d df", dc_fmt3(pr$wald), pr$wald_df))
            add(3, "Probability that the test rejects under parallel trends (size)", pr$size, if (k > 1) sprintf("with %d pre-treatment effects the chance of a false alarm is above the nominal %s%%", k, dc_fmt3(100 - o$ciLevel)) else "")
            if (is.finite(pr$power)) add(4, sprintf("Power against a trend of %s per period", dc_fmt3(pr$slope)), pr$power, sprintf("a violation growing %s per period would be detected with this probability", dc_fmt3(pr$slope)))
            add(5, "Slope detected with 50% power", pr$slope_for_power[1], "per event period; the violation at event time e is slope x (e + 1)")
            add(6, "Slope detected with 80% power", pr$slope_for_power[2], "")
            pt <- self$results$ptPost
            for (i in seq_along(pr$e_post)) {
                e <- pr$e_post[i]; rel <- e - pr$ref
                vals <- list(e = e, est = pr$beta_post[i], se = pr$se_post[i], b50 = pr$slope_for_power[1] * rel, b80 = pr$slope_for_power[2] * rel)
                if (is.finite(pr$power)) vals$bh <- pr$slope * rel
                pt$addRow(rowKey = i, values = vals)
            }
            pt$setNote("a", "If the bias in a post-treatment effect is as large as the estimated effect, a violation of parallel trends that the pre-test would usually miss could account for the whole effect.")
            self$results$ptPlot$setState(list(pr = pr))
        },

        .plotPt = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            pr <- st$pr
            es <- data.frame(e = c(pr$e_pre, pr$e_post), est = c(pr$beta_pre, pr$beta_post), se = c(pr$se_pre, pr$se_post))
            es <- rbind(es, data.frame(e = pr$ref, est = 0, se = 0))
            xs <- seq(min(es$e), max(es$e), length.out = 50)
            lines <- data.frame(e = rep(xs, 2 + is.finite(pr$power)),
                                trend = c(pr$slope_for_power[1] * (xs - pr$ref), pr$slope_for_power[2] * (xs - pr$ref), if (is.finite(pr$power)) pr$slope * (xs - pr$ref)),
                                what = rep(c("Detected with 50% power", "Detected with 80% power", if (is.finite(pr$power)) "Hypothesized trend"), each = length(xs)))
            g <- ggplot2::ggplot() + ggtheme +
                ggplot2::geom_hline(yintercept = 0, colour = "grey60") + ggplot2::geom_vline(xintercept = pr$ref + 0.5, linetype = "dotted", colour = "grey60") +
                ggplot2::geom_line(data = lines, ggplot2::aes(x = e, y = trend, colour = what, linetype = what), linewidth = 0.8) +
                ggplot2::geom_errorbar(data = es[es$se > 0, ], ggplot2::aes(x = e, ymin = est - pr$crit * se, ymax = est + pr$crit * se), width = 0.15, colour = "grey30") +
                ggplot2::geom_point(data = es, ggplot2::aes(x = e, y = est), size = 2.2, colour = "grey20") +
                ggplot2::labs(x = "Event time", y = "Estimated effect", colour = NULL, linetype = NULL, caption = "Lines: linear violations of parallel trends through the reference period that the pre-test rejects with the stated probability")
            print(g)
            TRUE
        },

        .plotBacon = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            print(dc_plot_bacon(s$detail, s$twfe, ggtheme, self$options$plotPalette))
            TRUE
        },

        .plotHonest = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            print(dc_plot_honest(s, ggtheme, self$options$plotPalette))
            TRUE
        },

        .notesHtml = function(o, facts) {
            n <- character()
            if (o$bacon) {
                n <- c(n, paste0("<b>Goodman-Bacon decomposition</b>: the static TWFE coefficient is a weighted average of all ",
                    "two-group, two-period comparisons. Comparisons of <i>later</i> against <i>earlier</i> treated units use ",
                    "already-treated units as controls and can bias the estimate when effects grow over time."))
                if (!is.null(facts$bacon))
                    n <- c(n, sprintf("Here %.0f%% of the TWFE weight is on later-vs-earlier-treated comparisons.",
                                      100 * facts$bacon$bad))
            }
            if (o$weights) {
                n <- c(n, paste0("<b>TWFE weights</b>: the TWFE coefficient equals a weighted sum of the effects of all treated ",
                    "observations. Negative weights mean some treated observations enter with the wrong sign, so the coefficient ",
                    "need not be a sensible average of effects."))
                if (!is.null(facts$weights))
                    n <- c(n, sprintf("Here %d of %d treated observations (%.1f%%) carry negative weight, summing to %.3f.",
                                      facts$weights$neg_n, facts$weights$n,
                                      100 * facts$weights$neg_n / facts$weights$n, facts$weights$neg_w))
            }
            if (o$pretest)
                n <- c(n, paste0("<b>Pre-trend tests</b> are joint Wald tests that all pre-treatment effects are zero. TWFE pre-trends ",
                    "can be contaminated by effect heterogeneity just like its post-treatment estimates. Not rejecting does not prove ",
                    "parallel trends: pre-tests can have low power and passing them can distort later inference (Roth 2022). ",
                    "The Callaway-Sant'Anna test uses all pre-treatment group-time effects regardless of the event window."))
            if (o$honest)
                n <- c(n, paste0("<b>HonestDiD</b> (Rambachan &amp; Roth 2023) reports confidence intervals that stay valid if parallel trends ",
                    "fail after treatment by up to Mbar times the largest violation seen before treatment. Mbar = 0 assumes parallel trends hold; ",
                    "the value of Mbar at which the interval first includes zero shows how robust the finding is."))
            if (o$ptpower)
                n <- c(n, paste0("<b>Pre-trend power</b> (Roth 2022): the usual pre-test passes if no pre-treatment effect is individually significant. ",
                    "This reports the probability that it would reject if parallel trends failed by a linear trend through the reference period (the bias at event time e is slope x (e + 1)), ",
                    "and the slopes that the test detects with 50 and 80 percent power. A pre-test with low power passes even when the trend is large enough to explain the whole effect; ",
                    "with several pre-treatment effects it also rejects too often when trends are parallel."))
            paste0("<p>", paste(n, collapse = "</p><p>"), "</p>")
        })
)
