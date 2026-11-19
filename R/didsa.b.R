
didsaClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "didsaClass",
    inherit = didsaBase,
    private = list(
        .run = function() dc_single_run(self, "sa"),
        .plot = function(image, ggtheme, theme, ...) dc_single_plot(image, ggtheme, self$options)
    )
)
