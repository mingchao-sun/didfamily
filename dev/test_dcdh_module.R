# The jamovi analysis "Treatment Switching On and Off" (R/dcdh.b.R) run headless on the sample data and on edge cases.
# Run from the module root, ideally in jamovi's own R:  Rscript dev/test_dcdh_module.R
root <- if (file.exists("R/core_dcdh.R")) "." else ".."
suppressMessages({library(jmvcore); library(R6); library(ggplot2)})
for (f in c("00jmv.R", "core_dcdh.R", "core_guide.R", "dcdh.h.R", "dcdh.b.R")) source(file.path(root, "R", f))
d <- read.csv(file.path(root, "dev/csv/sim_dcdh.csv"))
bad <- 0
run <- function(label, data, expect_error = FALSE, ...) {
    r <- tryCatch({ res <- dcdh(data = data, ...); invisible(capture.output(print(res))); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
    res <- NULL
    good <- if (expect_error) grepl("^ERROR", r) else r == "ok"
    cat(sprintf("%-46s %s\n", label, if (good) substr(r, 1, 80) else paste("UNEXPECTED:", r))); if (!good) bad <<- bad + 1
    invisible(r)
}
run("sample data, defaults", d, dep = "y", id = "id", time = "year", treat = "D")
run("guide, notes, gray plot", d, dep = "y", id = "id", time = "year", treat = "D", showGuide = TRUE, showHelp = TRUE, plotPalette = "gray")
run("5 effects, 3 placebos, 90%", d, dep = "y", id = "id", time = "year", treat = "D", effects = 5, placebo = 3, level = 90)
run("no placebos", d, dep = "y", id = "id", time = "year", treat = "D", placebo = 0)
d2 <- d; d2$D <- factor(d2$D)
run("nominal treatment", d2, dep = "y", id = "id", time = "year", treat = "D")
d3 <- d[-sample(nrow(d), 150), ]
run("unbalanced panel", d3, dep = "y", id = "id", time = "year", treat = "D")
d4 <- d; d4$D <- d4$D * 2
run("treatment not 0/1", d4, TRUE, dep = "y", id = "id", time = "year", treat = "D")
d5 <- d; d5$D <- 0
run("nobody switches", d5, TRUE, dep = "y", id = "id", time = "year", treat = "D")
d6 <- d; d6$year <- d6$year + 0.5
run("time not whole numbers", d6, TRUE, dep = "y", id = "id", time = "year", treat = "D")
run("missing treatment variable", d, dep = "y", id = "id", time = "year", treat = NULL)
d7 <- rbind(d, d[1:5, ])
run("duplicate group-period rows", d7, TRUE, dep = "y", id = "id", time = "year", treat = "D")
# numbers equal the engine
res <- dcdh(data = d, dep = "y", id = "id", time = "year", treat = "D", effects = 3, placebo = 2)
tab <- as.data.frame(res$effects)
me <- dc_dcdh(d$y, d$id, d$year, d$D, effects = 3, placebo = 2)
ok <- max(abs(tab$estimate - me$table$estimate)) < 1e-12 && max(abs(tab$se - me$table$se)) < 1e-12
cat(sprintf("%-46s %s\n", "table equals the engine", if (ok) "ok" else "FAIL")); if (!ok) bad <- bad + 1
cat(if (bad == 0) "\nAll checks passed.\n" else sprintf("\n%d check(s) failed.\n", bad))
if (bad > 0) quit(status = 1)
