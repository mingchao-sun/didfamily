
didcsClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didcsClass",
    inherit = didcsBase,
    private = list(
        .run = function() dc_single_run(self, "cs"),
        .plot = function(image, ggtheme, theme, ...) dc_single_plot(image, ggtheme, self$options),
        .plotGT = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(dc_plot_gt(image$state, ggtheme, self$options$plotPalette))
            TRUE
        }
    )
)
