# The pre-trend power option of the analysis "Bacon, Pre-trend and HonestDiD" (R/didiag.b.R) run headless on the sample data and on edge cases.
# Run from the module root, ideally in jamovi's own R:  Rscript dev/test_pretrend_module.R
root <- if (file.exists("R/core_pretrend.R")) "." else ".."
suppressMessages({library(jmvcore); library(R6); library(ggplot2)})
for (f in c("00jmv.R", "core_prepare.R", "core_estimators.R", "core_diag.R", "core_plot.R", "core_pretrend.R", "core_guide.R", "didiag.h.R", "didiag.b.R")) source(file.path(root, "R", f))
d <- read.csv(file.path(root, "dev/csv/sim_staggered.csv")); print(names(d))
bad <- 0
run <- function(label, data, expect_error = FALSE, ...) {
    out <- NULL
    r <- tryCatch({ res <- didiag(data = data, ...); out <<- res; invisible(capture.output(print(res))); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
    good <- if (expect_error) grepl("^ERROR", r) else r == "ok"
    cat(sprintf("%-52s %s\n", label, if (good) substr(r, 1, 80) else paste("UNEXPECTED:", r))); if (!good) bad <<- bad + 1
    invisible(out)
}
a <- list(dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", bacon = FALSE, weights = FALSE, pretest = FALSE, honest = FALSE)
ex <- function(...) c(modifyList(a, list(...)))
call <- function(label, data, ..., expect_error = FALSE) do.call(run, c(list(label, data, expect_error), ex(...)))
r1 <- call("Sun-Abraham event study, power for no slope", d, ptpower = TRUE)
r2 <- call("Callaway-Sant'Anna event study", d, ptpower = TRUE, ptBase = "cs")
r3 <- call("hypothesized slope 0.1, notes, guide", d, ptpower = TRUE, ptSlope = 0.1, showHelp = TRUE, showGuide = TRUE)
r4 <- call("event window with one pre-period", d, ptpower = TRUE, minE = -2)
r5 <- call("no pre-period beyond the reference (message in the table)", d, ptpower = TRUE, minE = -1)
r6 <- call("together with HonestDiD off and the pre-test on", d, ptpower = TRUE, pretest = TRUE)
if (!is.null(r1)) { tb <- r1$ptTable$asDF; cat("\n"); print(tb[, c("item", "value")]); stopifnot(nrow(tb) >= 5, all(is.finite(tb$value[!is.na(tb$value)]))) }
if (!is.null(r3)) { tb <- r3$ptTable$asDF; stopifnot(any(grepl("Power against a trend", tb$item))); pb <- r3$ptPost$asDF; stopifnot("bh" %in% names(pb)) }
cat(if (bad == 0) "\nAll checks passed.\n" else sprintf("\n%d check(s) failed.\n", bad))
if (bad > 0) quit(status = 1)
