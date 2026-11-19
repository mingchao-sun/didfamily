# Checks the variable types and descriptions stored in the .omv example datasets.
root <- if (file.exists("data/castle.omv")) "." else ".."
check <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) quit(status = 1) }
fields <- function(ds) {
    f <- file.path(root, "data", paste0(ds, ".omv"))
    meta <- jsonlite::fromJSON(unz(f, "metadata.json"), simplifyVector = FALSE)
    stats::setNames(meta$dataSet$fields, vapply(meta$dataSet$fields, function(x) x$name, ""))
}
expect <- list(
    sim_staggered = c(id = "ID", year = "Continuous", y = "Continuous", D = "Nominal", g = "Continuous"),
    sim_2x2 = c(id = "ID", year = "Continuous", y = "Continuous", treated = "Nominal", x = "Continuous", region = "Nominal"),
    sim_ddd = c(id = "ID", year = "Continuous", y = "Continuous", g = "Continuous", elig = "Nominal"),
    sim_dose = c(id = "ID", year = "Continuous", y = "Continuous", g = "Continuous", dose = "Continuous"),
    sim_dcdh = c(id = "ID", year = "Continuous", y = "Continuous", D = "Nominal"),
    sim_etwfe = c(id = "ID", year = "Continuous", y = "Continuous", g = "Continuous", x1 = "Continuous"),
    mpdta = c(countyreal = "ID", year = "Continuous", lemp = "Continuous", first.treat = "Continuous", treat = "Nominal", lpop = "Continuous"),
    castle = c(state = "ID", year = "Continuous", post = "Nominal", effyear = "Continuous", l_homicide = "Continuous", pop2000 = "Continuous"))
for (ds in names(expect)) {
    fl <- fields(ds)
    got <- vapply(names(expect[[ds]]), function(v) if (v %in% names(fl)) fl[[v]]$measureType else NA_character_, "")
    check(all(got == expect[[ds]]), sprintf("%s: variable types are as designed (%s)", ds, paste(names(got), got, sep = "=", collapse = ", ")))
    check(all(vapply(fl, function(x) nzchar(x$description), TRUE)), sprintf("%s: every variable has a description", ds))
}
# the datasets listed in the module metadata exist and nothing else (e.g. CSV) ships in data/
y <- yaml::read_yaml(file.path(root, "jamovi", "0000.yaml"))
paths <- vapply(y$datasets, function(x) x$path, "")
check(all(file.exists(file.path(root, "data", paths))) && all(grepl("[.]omv$", paths)), "all listed datasets exist and are .omv")
check(setequal(list.files(file.path(root, "data")), paths), "data/ contains exactly the listed datasets")
# the .omv data match the CSV data
for (ds in names(expect)) {
    csv <- read.csv(file.path(root, "dev", "csv", paste0(ds, ".csv")), check.names = FALSE)
    ref <- jmvReadWrite::read_omv(file.path(root, "data", paste0(ds, ".omv")))
    num <- intersect(names(csv), names(ref))
    tonum <- function(x) if (is.logical(x)) as.numeric(x) else suppressWarnings(as.numeric(as.character(x)))   # read_omv returns 0/1 as logical
    same <- vapply(num, function(v) isTRUE(all.equal(tonum(csv[[v]]), tonum(ref[[v]]), check.attributes = FALSE)) ||
                                      identical(as.character(csv[[v]]), as.character(ref[[v]])), TRUE)
    check(all(same) && nrow(csv) == nrow(ref), sprintf("%s: .omv values equal the CSV (%d rows)", ds, nrow(ref)))
}
cat("\nAll omv checks passed.\n")
