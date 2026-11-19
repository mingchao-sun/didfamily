# The six single-estimator analyses must give exactly the numbers of "Compare DiD Estimators".
# Run from the project root in a jamovi-like environment, with the module installed:
#   R_LIBS_USER=notthere R_LIBS="<module R>;<jamovi>/modules/jmv/R;<jamovi>/modules/base/R" Rscript dev/test_single.R
suppressMessages(library(didfamily))
ok <- 0L; bad <- 0L
chk <- function(what, cond) { if (isTRUE(cond)) ok <<- ok + 1L else { bad <<- bad + 1L; cat("FAIL:", what, "\n") } }

fun <- list(twfe = didtwfe, cs = didcs, sa = didsa, bjs = didbjs, gardner = didgardner, stacked = didstacked)
flag <- function(k) setNames(as.list(names(fun) == k), names(fun))

run_both <- function(d, args, label, skip = character()) {
    cmp <- do.call(didcompare, c(list(data = d), args))$summary$asDF
    for (k in setdiff(names(fun), skip)) {
        a <- args[setdiff(names(args), if (k == "bjs") "weights")]
        one <- do.call(fun[[k]], c(list(data = d), a, list(showCode = TRUE)))
        s <- one$summary$asDF
        row <- cmp[rownames(cmp) == paste0('"', k, '"'), ]
        chk(paste(label, k, "att"), isTRUE(all.equal(s$att, row$att, tolerance = 1e-8)))
        chk(paste(label, k, "se"),  isTRUE(all.equal(s$se, row$se, tolerance = 1e-8)))
        dy <- one$dynamic$asDF
        chk(paste(label, k, "dynamic table"), nrow(dy) >= 2 && nrow(dy) <= (args$maxE %||% 4) + (-(args$minE %||% -4)) + 1)
        chk(paste(label, k, "reference row"), any(dy$e == -1 & dy$est == 0))
        chk(paste(label, k, "code"), nchar(one$code$content) > 200)
    }
}
`%||%` <- function(a, b) if (is.null(a)) b else a

d <- read.csv("dev/csv/sim_staggered.csv")
run_both(d, list(dep = "y", id = "id", time = "year", treat = "D"), "sim/indicator")
run_both(d, list(dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", minE = -3, maxE = 5), "sim/cohort")

m <- read.csv("dev/csv/mpdta.csv")
run_both(m, list(dep = "lemp", id = "countyreal", time = "year", treat = "first.treat", treatType = "cohort", covs = "lpop"), "mpdta/covs")

ca <- read.csv("dev/csv/castle.csv")
run_both(ca, list(dep = "l_homicide", id = "state", time = "year", treat = "post", weights = "pop2000"), "castle/weights", skip = "bjs")   # the comparison skips BJS with weights; the BJS analysis has no weights option

# Callaway & Sant'Anna specifics
r <- didcs(data = d, dep = "y", id = "id", time = "year", treat = "D", csControl = "notyettreated",
           aggGroup = TRUE, aggCalendar = TRUE, showGT = TRUE)
c2 <- didcompare(data = d, dep = "y", id = "id", time = "year", treat = "D", csControl = "notyettreated",
                 twfe = FALSE, sa = FALSE, bjs = FALSE, gardner = FALSE, stacked = FALSE, aggGroup = TRUE, aggCalendar = TRUE)
chk("cs notyettreated att", isTRUE(all.equal(r$summary$asDF$att, c2$summary$asDF$att)))
chk("cs cohort table", isTRUE(all.equal(r$groupTab$asDF$est, c2$groupTab$asDF$est)))
chk("cs calendar table", nrow(r$calTab$asDF) == nrow(c2$calTab$asDF))
rb <- didcs(data = d, dep = "y", id = "id", time = "year", treat = "D", bootstrap = TRUE, bootReps = 199)
chk("cs bootstrap note", grepl("bootstrap", rb$summary$asDF$control))

# BJS ignores nothing silently: it has no weights option
chk("bjs has no weights option", !("weights" %in% names(formals(didbjs))))
# method notes and data guide: present in every analysis, content set whatever the selection
for (f in list(didtwfe, didcs, didsa, didbjs, didgardner, didstacked, didcompare, didiag, did22, dosedid, dddcompare, didsetup))
    chk("notes and guide options", all(c("showHelp", "showGuide") %in% names(formals(f))))
g <- didcs(data = d, dep = "y", id = "id", time = "year", treat = "D", showGuide = TRUE, showHelp = TRUE)
chk("guide text", grepl("Treatment variable", g$guide$content) && grepl("Sample data 1", g$guide$content) && grepl("This run", g$notes$content) && grepl("<h2>Data</h2>", g$notes$content))
chk("notes text", grepl("Callaway", g$notes$content))
# static method notes exist before any variable is chosen
for (k in c("didsetup", "did22", "didcompare", "didiag", "dosedid", "dddcompare", "didtwfe", "didcs", "didsa", "didbjs", "didgardner", "didstacked", "didetwfe")) {
    txt <- if (k %in% c("didtwfe", "didcs", "didsa", "didbjs", "didgardner", "didstacked", "didetwfe")) didfamily:::dc_single_notes(sub("^did", "", k), list()) else didfamily:::dc_method_html(k, list())
    chk(paste("static notes", k), nchar(txt) > 100)
    chk(paste("guide", k), nchar(didfamily:::dc_guide_html(k)) > 100)
}
cat(sprintf("single-estimator tests: %d passed, %d failed\n", ok, bad))
if (bad) quit(status = 1)
