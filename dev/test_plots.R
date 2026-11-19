# Renders every plot function with every palette / CI style and checks that it draws without
# errors or warnings. Optional 1st argument: a directory to save PNGs for visual inspection.
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_estimators", "core_did2s", "vendor_honestdid", "core_diag", "core_2x2", "core_dose", "core_panel", "core_plot"))
    source(file.path(root, "R", paste0(f, ".R")))
out_dir <- commandArgs(TRUE)[1]
check <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) quit(status = 1) }
draw <- function(p, name) {
    f <- if (is.na(out_dir)) tempfile(fileext = ".png") else file.path(out_dir, paste0(name, ".png"))
    warn <- character()
    ok <- tryCatch({ withCallingHandlers({ grDevices::png(f, width = 720, height = 450); print(p); grDevices::dev.off() },
                         warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }); TRUE },
                   error = function(e) { try(grDevices::dev.off(), silent = TRUE); cat("   error:", conditionMessage(e), "\n"); FALSE })
    if (length(warn)) cat("   warnings:", paste(unique(warn), collapse = " | "), "\n")
    check(ok && length(warn) == 0 && file.exists(f) && file.size(f) > 3000, paste("draws without errors or warnings:", name))
}

sim <- read.csv(file.path(root, "dev/csv/sim_staggered.csv"))
d <- dc_prepare(sim, "y", "id", "year", "D", "indicator")
res <- list(dc_twfe(d, -4, 4), dc_cs(d, -4, 4, extras = TRUE), dc_sa(d, -4, 4), dc_bjs(d, -4, 4), dc_gardner(d, -4, 4), dc_stacked(d, -4, 4))
z <- 1.96
long <- do.call(rbind, lapply(res, function(r) data.frame(method = r$method, e = r$dynamic$e, est = r$dynamic$est,
                                                         lo = r$dynamic$est - z * r$dynamic$se, hi = r$dynamic$est + z * r$dynamic$se)))
forest <- do.call(rbind, lapply(res, function(r) data.frame(method = r$method, est = r$att, lo = r$att - z * r$se, hi = r$att + z * r$se)))
for (pal in c("color", "gray")) {
    for (ci in c("bars", "ribbon"))
        draw(dc_plot_event(long, NULL, pal, ci, 0.6), paste0("event_", pal, "_", ci))
    draw(dc_plot_forest(forest, NULL, pal), paste0("forest_", pal))
    draw(dc_plot_gt(res[[2]]$extra$gt, NULL, pal), paste0("gt_", pal))
    bc <- dc_bacon(d); draw(dc_plot_bacon(bc$detail, bc$twfe, NULL, pal), paste0("bacon_", pal))
    draw(dc_plot_honest(data.frame(label = c("Original", "Mbar = 0", "Mbar = 1"), lb = c(.5, .5, .4), ub = c(.7, .7, .8), original = c(TRUE, FALSE, FALSE)), NULL, pal), paste0("honest_", pal))
    pc <- dc_panel_check(sim, "id", "year", "D", "indicator"); draw(dc_plot_status(pc$status, NULL, pal), paste0("status_", pal))
}
long2 <- long[long$method %in% c("TWFE", "Sun-Abraham"), ]
draw(dc_plot_event(long2, NULL, "gray", "bars", 0), "event_two_methods_nodge")
d22 <- dc_prepare22(read.csv(file.path(root, "dev/csv/sim_2x2.csv")), "y", "treated", "year", id = "id")
for (pal in c("color", "gray")) draw(dc_plot_22(dc_means22(d22), d22$labels, NULL, pal), paste0("twobytwo_", pal))
dd <- read.csv(file.path(root, "dev/csv/sim_dose.csv")); dp <- dc_prepare(dd, "y", "id", "year", "g", "cohort", dose = "dose")
rd <- dc_dose(dp, "nevertreated", 2, 0)
for (pal in c("color", "gray")) draw(dc_plot_dose(rd$curve, z, NULL, pal), paste0("dose_", pal))
# a real jamovi theme object must work too (its own colour scale is overridden by ours)
th <- jmvcore::theme_default()
draw(dc_plot_event(long, th, "gray", "ribbon"), "event_with_jamovi_theme")
cat("\nAll plot tests passed.\n")
