
didsetupClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didsetupClass",
    inherit = didsetupBase,
    private = list(

        .run = function() {
            self$results$guide$setContent(dc_guide_html("didsetup"))
            dc_hide_empty(self, c("id", "time", "treat"))
            self$results$notes$setContent(dc_method_html("didsetup", self$options))
            o <- self$options
            if (is.null(o$id) || is.null(o$time) || is.null(o$treat)) {
                self$results$info$setContent(paste0(
                    "<p>Select a unit ID, a time period and a treatment variable to check the panel before estimating. ",
                    "Problems are listed rather than stopping the analysis.</p>"))
                return()
            }
            r <- dc_panel_check(self$data, o$id, o$time, o$treat, o$treatType, dep = o$dep)

            tab <- self$results$checks
            label <- c(ok = "OK", note = "Note", warn = "Check", problem = "Problem")
            for (i in seq_len(nrow(r$checks)))
                tab$addRow(rowKey = i, values = list(
                    check = r$checks$check[i], status = label[[r$checks$status[i]]],
                    result = r$checks$result[i]))

            if (!is.null(r$cohorts)) {
                ct <- self$results$cohorts
                for (i in seq_len(nrow(r$cohorts))) {
                    vals <- list(cohort = r$cohorts$cohort[i], units = r$cohorts$units[i],
                                 share = r$cohorts$share[i])
                    if (!is.na(r$cohorts$pre[i])) { vals$pre <- r$cohorts$pre[i]; vals$post <- r$cohorts$post[i] }
                    ct$addRow(rowKey = i, values = vals)
                }
            }
            if (!is.null(r$status))
                self$results$plot$setState(list(status = r$status, sampled = r$n_sampled, total = r$n_units))
        },

        .plot = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            p <- dc_plot_status(s$status, ggtheme, self$options$plotPalette)
            if (s$sampled < s$total)
                p <- p + ggplot2::labs(caption = sprintf("Showing %d of %d units", s$sampled, s$total))
            print(p)
            TRUE
        })
)
