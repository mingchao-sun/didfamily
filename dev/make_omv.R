# Builds the .omv example datasets (variable types and descriptions preset) from the CSV files.
# The CSVs live in dev/csv/ (for scripts and tests); only the .omv files are in data/, which is what
# jamovi/0000.yaml lists and what ships with the module (CSVs in data/ would also be installed as R datasets).
# Needs the CRAN package jmvReadWrite (development only, not a module dependency):
#   Rscript dev/make_omv.R
suppressMessages(library(jmvReadWrite))
root <- if (file.exists("dev/csv/sim_staggered.csv")) "." else ".."
csv_dir <- file.path(root, "dev", "csv")
omv_dir <- file.path(root, "data")

# kind: "id" (ID), "cont" (continuous), "bin" (nominal 0/1, integer), "nom" (nominal text)
spec <- list(
    sim_staggered = list(
        id = c("id", "Unit identifier"),
        year = c("cont", "Period, 1 to 12"),
        y = c("cont", "Outcome"),
        D = c("bin", "Treatment indicator: 1 once the unit is treated (never turns off)"),
        g = c("cont", "First treated period of the unit, 0 = never treated")),
    sim_2x2 = list(
        id = c("id", "Unit identifier"),
        year = c("cont", "Period: 2010 = before, 2015 = after"),
        y = c("cont", "Outcome"),
        treated = c("bin", "Treated group: 1 = treated, 0 = comparison"),
        x = c("cont", "Baseline covariate; the outcome trend depends on it"),
        region = c("nom", "Region (north, south, east)")),
    sim_ddd = list(
        id = c("id", "Unit identifier"),
        year = c("cont", "Period, 1 to 8"),
        y = c("cont", "Outcome"),
        g = c("cont", "First treated period of the unit's group, 0 = never treated"),
        elig = c("bin", "Eligibility: 1 = eligible for the policy, 0 = not eligible")),
    sim_dose = list(
        id = c("id", "Unit identifier"),
        year = c("cont", "Period, 1 to 4"),
        y = c("cont", "Outcome"),
        g = c("cont", "First treated period of the unit, 0 = never treated"),
        dose = c("cont", "Dose received once treated (0 if untreated)")),
    sim_dcdh = list(
        id = c("id", "Group identifier"),
        year = c("cont", "Period, 2001 to 2010"),
        y = c("cont", "Outcome"),
        D = c("bin", "Treatment: 1 while the group is treated; it switches on in 218 groups and off again in about a third of them")),
    sim_etwfe = list(
        id = c("id", "Unit identifier"),
        year = c("cont", "Period, 1 to 8"),
        y = c("cont", "Count outcome"),
        g = c("cont", "First treated period of the unit, 0 = never treated"),
        x1 = c("cont", "Baseline covariate: it moves the level of the outcome and the size of the effect")),
    mpdta = list(
        countyreal = c("id", "County identifier"),
        year = c("cont", "Year, 2003 to 2007"),
        lemp = c("cont", "Log of county teen employment"),
        first.treat = c("cont", "Year of the first minimum-wage increase, 0 = never treated"),
        treat = c("bin", "1 = county is ever treated"),
        lpop = c("cont", "Log of county population")),
    castle = list(
        state = c("id", "US state"),
        year = c("cont", "Year, 2000 to 2010"),
        post = c("bin", "Treatment indicator: 1 once the state's castle doctrine law is in force"),
        effyear = c("cont", "Year the law took effect, 0 = never treated"),
        l_homicide = c("cont", "Log homicide rate"),
        l_robbery = c("cont", "Log robbery rate"),
        l_assault = c("cont", "Log assault rate"),
        l_burglary = c("cont", "Log burglary rate"),
        l_larceny = c("cont", "Log larceny rate"),
        l_motor = c("cont", "Log motor vehicle theft rate"),
        population = c("cont", "State population"),
        income = c("cont", "Median household income"),
        unemployrt = c("cont", "Unemployment rate, percent"),
        poverty = c("cont", "Poverty rate, percent"),
        police = c("cont", "Police officers per 100k residents"),
        pop2000 = c("cont", "State population in 2000: a time-constant sampling weight")))

for (ds in names(spec)) {
    d <- read.csv(file.path(csv_dir, paste0(ds, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(setequal(names(d), names(spec[[ds]])))
    d <- d[, names(spec[[ds]]), drop = FALSE]
    for (v in names(d)) {
        kind <- spec[[ds]][[v]][1]; desc <- spec[[ds]][[v]][2]
        if (kind == "id") {
            if (is.numeric(d[[v]])) d[[v]] <- as.integer(d[[v]])
            attr(d[[v]], "jmv-id") <- TRUE
        } else if (kind == "bin") {
            d[[v]] <- factor(as.integer(d[[v]]), levels = c(0L, 1L))
        } else if (kind == "nom") {
            d[[v]] <- factor(d[[v]])
        } else if (kind == "cont" && all(d[[v]] == round(d[[v]]), na.rm = TRUE)) {
            d[[v]] <- as.integer(d[[v]])
        }
        attr(d[[v]], "jmv-desc") <- desc
    }
    out <- file.path(omv_dir, paste0(ds, ".omv"))
    write_omv(d, out, frcWrt = TRUE)
    cat(sprintf("%-14s %6d rows  %2d vars -> %s\n", ds, nrow(d), ncol(d), basename(out)))
}
