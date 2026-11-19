
didtwfeClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didtwfeClass",
    inherit = didtwfeBase,
    private = list(
        .run = function() dc_single_run(self, "twfe"),
        .plot = function(image, ggtheme, theme, ...) dc_single_plot(image, ggtheme, self$options)
    )
)
