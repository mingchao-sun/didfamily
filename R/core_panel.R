# Core layer: panel diagnostics. Unlike dc_prepare(), nothing here stops on a problem: every
# issue is reported as a row of the checks table so users can see what to fix.

# returns list(checks, cohorts, status) where
#   checks  : data.frame(check, result, status) with status in "ok", "note", "warn", "problem"
#   cohorts : data.frame(cohort, units, share, pre, post)
#   status  : data.frame(unit, t, status) for the treatment-status plot (a sample of units)
dc_panel_check <- function(data, id, time, treat, treat_type = c("indicator", "cohort"), dep = NULL,
                           max_units = 150) {
    treat_type <- match.arg(treat_type)
    checks <- data.frame(check = character(), result = character(), status = character(), stringsAsFactors = FALSE)
    add <- function(check, result, status = "ok")
        checks[nrow(checks) + 1, ] <<- c(check, result, status)

    n <- nrow(data)
    idv <- as.character(data[[id]]); tv <- dc_num(data[[time]]); rv <- dc_num(data[[treat]])
    yv <- if (!is.null(dep)) dc_num(data[[dep]]) else NULL

    # --- missing values ---
    miss <- c(id = sum(is.na(data[[id]])), time = sum(is.na(tv)), treatment = sum(is.na(rv)),
              if (!is.null(yv)) c(outcome = sum(is.na(yv))))
    miss_tv <- sum(is.na(tv) & !is.na(data[[time]]))
    if (treat_type == "cohort") miss["treatment"] <- 0      # missing means never treated
    if (any(miss > 0))
        add("Missing values", paste(sprintf("%s: %d", names(miss)[miss > 0], miss[miss > 0]), collapse = ", "), "warn")
    else add("Missing values", "none")
    keep <- !is.na(idv) & !is.na(tv)
    if (treat_type == "indicator") keep <- keep & !is.na(rv)
    if (!is.null(yv)) keep <- keep & !is.na(yv)
    if (miss_tv > 0) add("Time variable", "some values are not numeric", "problem")
    idv <- idv[keep]; tv <- tv[keep]; rv <- rv[keep]
    if (!length(idv)) {
        add("Usable rows", "none", "problem")
        return(list(checks = checks, cohorts = NULL, status = NULL))
    }

    periods <- sort(unique(tv)); units <- unique(idv)
    add("Size", sprintf("%d usable rows, %d units, %d periods (%g to %g)", length(idv), length(units), length(periods), min(periods), max(periods)))
    if (any(tv != round(tv))) add("Time values", "not integer-valued; use a period number or year", "problem")
    else {
        dd <- diff(periods)
        gaps <- setdiff(seq(min(periods), max(periods)), periods)
        if (length(dd) && all(dd == dd[1]) && dd[1] > 1)
            add("Time values", sprintf("equally spaced, %g apart; recoded to 1, 2, 3, ... for estimation (event time counts periods)", dd[1]), "note")
        else if (length(gaps))
            add("Time values", sprintf("not equally spaced (missing from the sequence: %s); event time is in time units", paste(head(gaps, 6), collapse = ", ")), "warn")
        else add("Time values", "consecutive integers")
        if (min(periods) <= 0)
            add("First period", sprintf("%g; recoded to start at 1 for estimation (0 is reserved for never treated)", min(periods)), "note")
    }

    dup <- sum(duplicated(data.frame(idv, tv)))
    if (dup) add("Duplicates", sprintf("%d repeated unit-period rows", dup), "problem") else add("Duplicates", "none")

    cnt <- table(idv)
    balanced <- all(cnt == length(periods)) && !dup
    if (balanced) add("Balance", "balanced panel (every unit observed in every period)")
    else add("Balance", sprintf("unbalanced: %d of %d units miss at least one period", sum(cnt < length(periods)), length(units)), "note")

    # --- treatment ---
    ord <- order(idv, tv); idv <- idv[ord]; tv <- tv[ord]; rv <- rv[ord]
    if (treat_type == "indicator") {
        if (!all(rv %in% c(0, 1))) {
            add("Treatment", "indicator must contain only 0 and 1", "problem")
            return(list(checks = checks, cohorts = NULL, status = NULL))
        }
        first <- tapply(ifelse(rv == 1, tv, Inf), idv, min)
        g <- as.numeric(first[idv]); g[is.infinite(g)] <- 0
        off <- tapply(seq_along(rv), idv, function(i) any(diff(rv[i]) < 0))
        n_off <- sum(off, na.rm = TRUE)
        if (n_off) add("Absorbing treatment", sprintf("%d units switch back to untreated; the estimators here need treatment that never turns off", n_off), "problem")
        else add("Absorbing treatment", "treatment never turns off")
        if (!balanced) add("Cohort definition", "first treated period is inferred from observed rows; a missing first treated row would shift a unit's cohort", "note")
        stat <- ifelse(rv == 1, "Treated", "Untreated")
    } else {
        g <- rv; g[is.na(g) | is.infinite(g)] <- 0
        varies <- tapply(g, idv, function(v) length(unique(v)) > 1)
        if (any(varies)) add("Cohort variable", sprintf("varies within %d units; it must be constant within a unit", sum(varies)), "problem")
        else add("Cohort variable", "constant within units")
        stat <- ifelse(g > 0 & tv >= g, "Treated", "Untreated")
    }
    tmin <- min(tv); tmax <- max(tv)
    unit_g <- tapply(g, idv, function(v) v[1])
    n_never <- sum(unit_g == 0); n_always <- sum(unit_g > 0 & unit_g <= tmin); n_late <- sum(unit_g > tmax)
    treated_units <- unit_g[unit_g > tmin & unit_g <= tmax]
    add("Never-treated units", sprintf("%d (%.0f%%)", n_never, 100 * n_never / length(unit_g)), if (n_never == 0) "note" else "ok")
    if (n_always) add("Treated from the first period", sprintf("%d units have no pre-period and are dropped by the estimators", n_always), "warn")
    if (n_late) add("Treated after the sample", sprintf("%d units are first treated after the last period and count as never treated", n_late), "note")
    cohorts_n <- table(treated_units)
    if (!length(cohorts_n)) add("Treatment timing", "no unit is treated inside the sample", "problem")
    else {
        add("Treatment timing", if (length(cohorts_n) == 1) "single treatment period: a 2x2 design is enough (see Classic DiD)"
            else sprintf("staggered adoption in %d cohorts", length(cohorts_n)), "ok")
        if (min(cohorts_n) < 10) add("Small cohorts", sprintf("smallest cohort has %d units; its estimates will be noisy", min(cohorts_n)), "warn")
    }

    # --- cohort table ---
    coh <- sort(unique(unit_g[unit_g > 0]))
    cohorts <- data.frame(cohort = c(if (n_never) "Never treated", as.character(coh)),
                          units = c(if (n_never) n_never, as.integer(table(factor(unit_g[unit_g > 0], levels = coh)))),
                          stringsAsFactors = FALSE)
    cohorts$share <- cohorts$units / length(unit_g)
    cohorts$pre <- c(if (n_never) NA_real_, pmax(coh - tmin, 0))
    cohorts$post <- c(if (n_never) NA_real_, pmax(tmax - coh + 1, 0))

    # --- treatment-status plot data: a sample of units sorted by cohort ---
    ukey <- data.frame(unit = names(unit_g), g = as.numeric(unit_g), stringsAsFactors = FALSE)
    ukey$sortg <- ifelse(ukey$g == 0, Inf, ukey$g)
    ukey <- ukey[order(ukey$sortg, ukey$unit), ]
    if (nrow(ukey) > max_units)
        ukey <- ukey[unique(round(seq(1, nrow(ukey), length.out = max_units))), ]
    ukey$row <- seq_len(nrow(ukey))
    grid <- expand.grid(unit = ukey$unit, t = periods, stringsAsFactors = FALSE)
    obs <- data.frame(unit = idv, t = tv, status = stat, stringsAsFactors = FALSE)
    grid <- merge(grid, obs, by = c("unit", "t"), all.x = TRUE)
    grid$status[is.na(grid$status)] <- "Not observed"
    grid$row <- ukey$row[match(grid$unit, ukey$unit)]
    list(checks = checks, cohorts = cohorts, status = grid[, c("row", "t", "status")],
         n_sampled = nrow(ukey), n_units = length(unit_g))
}
