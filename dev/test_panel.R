# Tests for panel diagnostics (plain R is enough: no estimators involved).
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_panel")) source(file.path(root, "R", paste0(f, ".R")))
check <- function(cond, msg) {
    cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
    if (!isTRUE(cond)) quit(status = 1)
}
st <- function(x, what) x$checks$status[x$checks$check == what]
rs <- function(x, what) x$checks$result[x$checks$check == what]

sim <- read.csv(file.path(root, "dev/csv/sim_staggered.csv"))
r <- dc_panel_check(sim, "id", "year", "D", "indicator", dep = "y")
check(grepl("600 units, 12 periods", rs(r, "Size")), "size summary")
check(st(r, "Balance") == "ok", "balanced panel recognised")
check(st(r, "Absorbing treatment") == "ok", "absorbing treatment recognised")
check(grepl("3 cohorts", rs(r, "Treatment timing")), "staggered adoption in 3 cohorts")
check(setequal(r$cohorts$cohort, c("Never treated", "5", "7", "9")) && abs(sum(r$cohorts$share) - 1) < 1e-12, "cohort table covers all units")
check(r$cohorts$units[r$cohorts$cohort == "5"] == 170, "cohort sizes match")
check(r$n_sampled == 150 && r$n_units == 600, "plot samples 150 of 600 units")
check(nrow(r$status) == 150 * 12 && all(r$status$status %in% c("Treated", "Untreated", "Not observed")), "status grid is units x periods")
first_rows <- r$status[r$status$row == 1, ]
check(all(first_rows$status[first_rows$t < 5] == "Untreated") && all(first_rows$status[first_rows$t >= 5] == "Treated"),
      "first plotted unit is a cohort-5 unit with the right status path")

rc <- dc_panel_check(sim, "id", "year", "g", "cohort", dep = "y")
check(rc$cohorts$units[rc$cohorts$cohort == "5"] == 170 && st(rc, "Cohort variable") == "ok", "cohort input gives the same cohorts")

# problems are reported, not thrown
bad <- sim; bad$D[bad$id == 3 & bad$year == 12] <- 0
bad$D[bad$id == which(sim$g[sim$year == 1] == 5)[1] & bad$year == 12] <- 0
rb <- dc_panel_check(bad, "id", "year", "D", "indicator")
check("problem" %in% st(rb, "Absorbing treatment") || any(grepl("switch back", rb$checks$result)), "non-absorbing treatment reported")
dupd <- rbind(sim, sim[1:3, ])
check(st(dc_panel_check(dupd, "id", "year", "D", "indicator"), "Duplicates") == "problem", "duplicate rows reported")
ub <- sim[-c(10, 500, 2000), ]
check(st(dc_panel_check(ub, "id", "year", "D", "indicator"), "Balance") == "note", "unbalanced panel reported")
miss <- sim; miss$y[1:5] <- NA; miss$year[7] <- NA
check(st(dc_panel_check(miss, "id", "year", "D", "indicator", dep = "y"), "Missing values") == "warn", "missing values reported")
bin <- sim; bin$D[1] <- 3
rbad <- dc_panel_check(bin, "id", "year", "D", "indicator")
check(st(rbad, "Treatment") == "problem" && is.null(rbad$cohorts), "invalid indicator reported and the function still returns")
one <- sim; one$g[] <- ifelse(one$g > 0, 5, 0); one$D <- as.integer(one$g > 0 & one$year >= 5)
check(grepl("single treatment period", rs(dc_panel_check(one, "id", "year", "D", "indicator"), "Treatment timing")), "single treatment period suggests a 2x2")
mp <- read.csv(file.path(root, "dev/csv/mpdta.csv"))
rm <- dc_panel_check(mp, "countyreal", "year", "first.treat", "cohort", dep = "lemp")
check(grepl("500 units, 5 periods", rs(rm, "Size")) && rm$n_sampled == 150, "mpdta summary")
cat("\nAll panel-diagnostics tests passed.\n")
