
didgardnerClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didgardnerClass",
    inherit = didgardnerBase,
    private = list(
        .run = function() dc_single_run(self, "gardner"),
        .plot = function(image, ggtheme, theme, ...) dc_single_plot(image, ggtheme, self$options)
    )
)
