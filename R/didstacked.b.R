
didstackedClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didstackedClass",
    inherit = didstackedBase,
    private = list(
        .run = function() dc_single_run(self, "stacked"),
        .plot = function(image, ggtheme, theme, ...) dc_single_plot(image, ggtheme, self$options)
    )
)
