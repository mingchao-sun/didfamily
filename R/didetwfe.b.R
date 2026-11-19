
didetwfeClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didetwfeClass",
    inherit = didetwfeBase,
    private = list(
        .run = function() dc_single_run(self, "etwfe"),
        .plot = function(image, ggtheme, theme, ...) dc_single_plot(image, ggtheme, self$options)
    )
)
