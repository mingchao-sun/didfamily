
dcdhClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "dcdhClass",
    inherit = dcdhBase,
    private = list(

        .run = function() {
            o <- self$options
            self$results$guide$setContent(dc_guide_html("dcdh"))
            self$results$notes$setContent(dc_method_html("dcdh", o))
            if (is.null(o$dep) || is.null(o$id) || is.null(o$time) || is.null(o$treat)) {
                self$results$info$setContent(paste0(
                    "<p>Select an outcome, a group ID, a time period and a binary treatment. The treatment may switch on and off; ",
                    "groups that never change it serve as comparisons until they do.</p>"))
                for (nm in c("effects", "tests", "plot")) self$results[[nm]]$setVisible(FALSE)
                return()
            }
            dat <- self$data
            y <- suppressWarnings(as.numeric(dat[[o$dep]]))
            tv <- dat[[o$treat]]
            d <- if (is.factor(tv)) {
                lv <- levels(droplevels(tv))
                if (length(lv) != 2) jmvcore::reject("The treatment must have exactly two values (the second one is treated).", code = "")
                as.numeric(droplevels(tv) == lv[2])
            } else suppressWarnings(as.numeric(as.character(tv)))
            tm <- suppressWarnings(as.numeric(as.character(dat[[o$time]])))
            g <- as.character(dat[[o$id]])
            if (any(is.finite(tm) & tm != round(tm))) jmvcore::reject("The time period must be whole numbers (a year or a period number).", code = "")
            n_all <- length(y)
            ok <- is.finite(y) & is.finite(d) & is.finite(tm) & !is.na(g)
            if (sum(ok) < 20) jmvcore::reject("Fewer than 20 complete observations.", code = "")
            if (any(!(d[ok] %in% c(0, 1)))) jmvcore::reject("The treatment must be binary (0/1).", code = "")
            res <- tryCatch(dc_dcdh(y[ok], g[ok], tm[ok], d[ok], effects = o$effects, placebo = o$placebo, level = o$level),
                            error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
            tab <- res$table
            ds <- res$design

            # ---- panel summary
            times <- sort(unique(tm[ok]))
            nin <- sum(is.finite(ds$F) & ds$d1 == 0); nout <- sum(is.finite(ds$F) & ds$d1 == 1)
            miss <- n_all - sum(ok)
            msg <- sprintf("<p>%d groups and %d periods. %d groups change their treatment (%d switch on from 0, %d switch off from 1); %d never change it.",
                           res$G, res$n_periods, res$n_switchers, nin, nout, res$G - res$n_switchers)
            if (miss > 0) msg <- paste0(msg, sprintf(" %d row(s) with missing values were dropped.", miss))
            if (any(diff(times) != 1)) msg <- paste0(msg, " Some periods are missing or unevenly spaced; the periods present are treated as consecutive, so check that this is what you intend.")
            if (res$G * res$n_periods > sum(ok)) msg <- paste0(msg, " The panel is unbalanced; groups are compared only in periods where they are observed.")
            if (any(tab$n_switchers == 0)) msg <- paste0(msg, " Some effects or placebos cannot be estimated: no group has enough periods or comparison groups.")
            self$results$info$setContent(paste0(msg, "</p>"))

            # ---- table
            et <- self$results$effects
            for (i in seq_len(nrow(tab))) {
                p <- if (is.finite(tab$se[i]) && tab$se[i] > 0) 2 * stats::pnorm(-abs(tab$estimate[i] / tab$se[i])) else NA_real_
                et$addRow(rowKey = i, values = list(
                    label = tab$label[i], lag = if (tab$kind[i] == "effect") sprintf("%d", tab$l[i] - 1) else sprintf("-%d", tab$l[i]),
                    estimate = tab$estimate[i], se = tab$se[i], lo = tab$lo[i], hi = tab$hi[i], p = p, n = tab$n_switchers[i]))
            }
            et$setNote("a", paste0("Each effect compares, for groups that change their treatment, the change in the outcome from the period before the first change ",
                                   "with the change of groups that have the same starting treatment and have not changed it yet. Effect 1 is the change in the period of the first switch (0 periods after it); ",
                                   "placebos look the same number of periods before the switch. Effects of groups switching off are signed so that they estimate the effect of the treatment."))

            tt <- self$results$tests
            jp <- dc_dcdh_joint(res, "placebo"); je <- dc_dcdh_joint(res, "effect")
            k <- 0
            if (!is.null(jp)) { k <- k + 1; tt$addRow(rowKey = k, values = list(test = "All placebos are zero", chi = jp$chi2, df = jp$df, p = jp$p)) }
            if (!is.null(je)) { k <- k + 1; tt$addRow(rowKey = k, values = list(test = "All effects are zero", chi = je$chi2, df = je$df, p = je$p)) }
            if (k == 0) tt$setVisible(FALSE)

            if (isTRUE(o$plot)) {
                st <- tab[is.finite(tab$estimate), , drop = FALSE]
                st$x <- ifelse(st$kind == "effect", st$l - 1, -st$l)
                self$results$plot$setState(st)
            }
        },

        .plot = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st) || nrow(st) == 0) return(FALSE)
            pal <- if (identical(self$options$plotPalette, "gray")) c(effect = "grey15", placebo = "grey55") else c(effect = "#0072B2", placebo = "#D55E00")
            st$kind <- factor(st$kind, levels = c("placebo", "effect"), labels = c("Placebo", "Effect"))
            names(pal) <- c("Effect", "Placebo")
            g <- ggplot2::ggplot(st, ggplot2::aes(x = x, y = estimate, colour = kind, shape = kind)) + ggtheme +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
                ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.15, na.rm = TRUE) +
                ggplot2::geom_point(size = 2.6) +
                ggplot2::scale_colour_manual(values = pal, name = NULL) + ggplot2::scale_shape_manual(values = c(Effect = 16, Placebo = 17), name = NULL) +
                ggplot2::scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]))) +
                ggplot2::labs(x = "Periods from the first change in treatment", y = "Estimate")
            print(g)
            TRUE
        })
)
