# Emulates jamovi's formula sandbox in plain R: stats::as.formula is replaced by
# jmvcore::asFormula, which rejects function calls outside a whitelist in formulas.
# Source this (after setting the jamovi-like R_LIBS) before running estimators.
dc_enable_sandbox <- function() {
    options(.__jmv.original_as.formula = stats::as.formula)
    patched <- function(object, env = parent.frame(), ...) jmvcore::asFormula(object, env)
    for (e in list(asNamespace("stats"), as.environment("package:stats"))) {
        unlockBinding("as.formula", e)
        assign("as.formula", patched, envir = e)
        lockBinding("as.formula", e)
    }
    invisible(TRUE)
}
