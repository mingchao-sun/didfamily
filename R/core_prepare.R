# Core layer: data standardisation. No jamovi dependencies in this file.
#
# Output columns (data frame `data`):
#   dc_y   outcome
#   dc_id  integer unit id
#   dc_t   integer-valued period
#   dc_g   first treated period (0 = never treated within the sample)
#   dc_d   treatment indicator (1 once treated)
#   dc_rel event time t - g (-1000 for never treated)

dc_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# count and first example of values that are present but not numbers after dc_num(); NULL if none
dc_first_bad <- function(raw, num) {
    r <- as.character(raw); bad <- !is.na(raw) & nzchar(trimws(r)) & is.na(num)
    if (any(bad)) c(sum(bad), r[which(bad)[1]]) else NULL
}

# covariates: numeric as is, factors/text as 0/1 dummies for every level but the first.
# Returns the new columns (dc_x1, dc_x2, ...) and a map code -> label / source variable / level.
dc_build_covs <- function(data, covs) {
    cols <- list()
    map <- data.frame(code = character(), label = character(), var = character(),
                      level = character(), stringsAsFactors = FALSE)
    for (v in covs) {
        x <- data[[v]]
        if (is.factor(x) || is.character(x)) {
            f <- factor(x)
            for (lev in levels(f)[-1]) {
                code <- sprintf("dc_x%d", nrow(map) + 1)
                cols[[code]] <- ifelse(is.na(f), NA_real_, as.numeric(f == lev))
                map[nrow(map) + 1, ] <- c(code, sprintf("%s = %s", v, lev), v, lev)
            }
        } else {
            code <- sprintf("dc_x%d", nrow(map) + 1)
            cols[[code]] <- dc_num(x)
            map[nrow(map) + 1, ] <- c(code, v, v, NA)
        }
    }
    list(cols = cols, map = map)
}

dc_prepare <- function(data, y, id, time, treat, treat_type = c("indicator", "cohort"), covs = NULL, part = NULL, cluster = NULL, dose = NULL, weights = NULL) {
    treat_type <- match.arg(treat_type)

    df <- data.frame(
        dc_y   = dc_num(data[[y]]),
        dc_id  = as.character(data[[id]]),
        dc_t   = dc_num(data[[time]]),
        dc_raw = dc_num(data[[treat]]),
        stringsAsFactors = FALSE)

    first_bad <- dc_first_bad
    b <- first_bad(data[[time]], df$dc_t)
    if (!is.null(b))
        stop(sprintf(paste0("The time variable must be numeric (a year or a period number): %s values such as '%s' are not numbers. ",
                            "Dates and text are not supported; convert them to a period number first."), b[1], b[2]))
    b <- first_bad(data[[y]], df$dc_y)
    if (!is.null(b))
        stop(sprintf("The outcome variable must be numeric: %s values such as '%s' are not numbers.", b[1], b[2]))
    b <- first_bad(data[[treat]], df$dc_raw)
    if (!is.null(b))
        stop(sprintf(if (treat_type == "indicator")
                         "The treatment indicator must be numeric 0/1: %s values such as '%s' are not numbers. Recode it to 0 and 1."
                     else "The first treated period must be a number (0 or empty for never treated): %s values such as '%s' are not numbers.",
                     b[1], b[2]))
    if (!is.null(part))
        df$dc_p <- dc_num(data[[part]])
    if (!is.null(cluster))
        df$dc_clraw <- as.character(data[[cluster]])
    if (!is.null(dose))
        df$dc_dose <- dc_num(data[[dose]])
    if (!is.null(weights)) {
        df$dc_sw <- dc_num(data[[weights]])
        b <- first_bad(data[[weights]], df$dc_sw)
        if (!is.null(b))
            stop(sprintf("The weight variable must be numeric: %s values such as '%s' are not numbers.", b[1], b[2]))
    }

    cb <- dc_build_covs(data, covs)
    for (cd in names(cb$cols)) df[[cd]] <- cb$cols[[cd]]
    cov_map <- cb$map

    df$dc_y[is.infinite(df$dc_y)] <- NA      # infinite outcomes are treated as missing
    n_input <- nrow(df)
    keep <- !is.na(df$dc_y) & !is.na(df$dc_id) & !is.na(df$dc_t)
    if (!is.null(part))
        keep <- keep & !is.na(df$dc_p)
    if (!is.null(cluster))
        keep <- keep & !is.na(df$dc_clraw)
    if (!is.null(weights))
        keep <- keep & !is.na(df$dc_sw)
    if (!is.null(dose))
        keep <- keep & (!is.na(df$dc_dose) | (treat_type == "cohort" & (is.na(df$dc_raw) | df$dc_raw == 0)))
    if (nrow(cov_map))
        keep <- keep & stats::complete.cases(df[, cov_map$code, drop = FALSE])
    if (treat_type == "indicator")
        keep <- keep & !is.na(df$dc_raw)
    df <- df[keep, , drop = FALSE]
    n_missing <- n_input - nrow(df)
    if (nrow(df) == 0)
        stop("No complete observations.")

    if (any(df$dc_t != round(df$dc_t)))
        stop("Time variable must be integer-valued (for example a year or period number).")

    # Periods are recoded to 1, 2, 3, ... when they start at zero or below (0 is the code for
    # "never treated") or are equally spaced by more than one (survey waves): estimators rely on
    # the period before treatment, g - 1, existing. Results are shown with the original values.
    periods <- sort(unique(df$dc_t))
    dd <- diff(periods)
    equal_step <- length(dd) > 0 && all(dd == dd[1])
    step <- if (equal_step) dd[1] else 1
    t0 <- min(periods)
    shifted <- t0 <= 0 || (equal_step && step > 1)
    time_map <- list(t0 = t0, step = step, shifted = shifted, irregular = length(dd) > 0 && !equal_step)
    if (shifted)
        df$dc_t <- (df$dc_t - t0) / step + 1
    df$dc_id <- match(df$dc_id, unique(df$dc_id))
    if (is.null(cluster)) {
        df$dc_cl <- df$dc_id
    } else {
        df$dc_cl <- match(df$dc_clraw, unique(df$dc_clraw))
        df$dc_clraw <- NULL
        if (any(stats::ave(df$dc_cl, df$dc_id, FUN = function(v) length(unique(v))) > 1))
            stop("The cluster variable must be constant within each unit (cluster at the unit level or above).")
    }

    if (anyDuplicated(df[, c("dc_id", "dc_t")]))
        stop("Some unit-period combinations appear more than once.")
    if (is.null(weights)) {
        df$dc_sw <- 1
    } else {
        if (any(df$dc_sw <= 0))
            stop("Weights must be positive.")
        if (any(stats::ave(df$dc_sw, df$dc_id, FUN = function(v) length(unique(v))) > 1))
            stop("The weight variable must be constant within each unit (for example a baseline population); time-varying weights are not supported.")
    }

    df <- df[order(df$dc_id, df$dc_t), , drop = FALSE]
    if (!is.null(part)) {
        if (!all(df$dc_p %in% c(0, 1)))
            stop("The eligibility (partition) variable must contain only 0 and 1.")
        if (any(stats::ave(df$dc_p, df$dc_id, FUN = function(v) length(unique(v))) > 1))
            stop("The eligibility (partition) variable must be constant within each unit.")
    }
    # covariates constant within every unit are absorbed by unit fixed effects
    cov_map$varying <- vapply(cov_map$code, function(cd)
        any(stats::ave(df[[cd]], df$dc_id, FUN = function(v) length(unique(v))) > 1), TRUE)
    tmin <- min(df$dc_t)
    tmax <- max(df$dc_t)

    if (treat_type == "indicator") {
        if (!all(df$dc_raw %in% c(0, 1)))
            stop("Treatment indicator must contain only 0 and 1.")
        first <- stats::ave(ifelse(df$dc_raw == 1, df$dc_t, Inf), df$dc_id, FUN = min)
        g <- ifelse(is.infinite(first), 0, first)
        if (any(g > 0 & df$dc_t >= g & df$dc_raw == 0))
            stop("Treatment is not absorbing: some units switch back to untreated. ",
                 "The estimators in this module need staggered adoption (treatment never turns off).")
    } else {
        g <- df$dc_raw
        g[is.na(g) | is.infinite(g)] <- 0
        if (any(stats::ave(g, df$dc_id, FUN = function(x) length(unique(x))) > 1))
            stop("The first-treated-period variable must be constant within each unit.")
        if (shifted) {
            pos <- g > 0
            g[pos] <- (g[pos] - t0) / step + 1
            if (any(pos & abs(g - round(g)) > 1e-8))
                stop(sprintf("Some first treated periods are not observed periods (periods start at %g and are %g apart).", t0, step))
            g[pos] <- pmax(g[pos], 1)      # treated before the sample starts: no pre-period, dropped below
        }
    }

    # treated only after the sample ends: indistinguishable from never treated
    g[g > tmax] <- 0
    # treated in the first observed period or before: no pre-period, drop
    always <- g > 0 & g <= tmin
    n_always <- length(unique(df$dc_id[always]))
    df$dc_g <- g
    if (!is.null(dose)) {
        df$dc_dose[is.na(df$dc_dose)] <- 0
        if (any(stats::ave(df$dc_dose, df$dc_id, FUN = function(v) length(unique(v))) > 1))
            stop("The dose variable must be constant within each unit (the dose it receives once treated).")
        if (any(df$dc_g > 0 & df$dc_dose <= 0))
            stop("Treated units need a positive dose.")
        df$dc_dose[df$dc_g == 0] <- 0
    }
    df <- df[!always, , drop = FALSE]
    if (nrow(df) == 0)
        stop("All units are treated in the first period.")
    if (!any(df$dc_g > 0))
        stop("No unit is treated inside the sample period (every unit is never treated), so there is nothing to estimate.")
    if (stats::var(df$dc_y) == 0)
        stop("The outcome variable is constant, so there is nothing to estimate.")

    df$dc_d   <- as.integer(df$dc_g > 0 & df$dc_t >= df$dc_g)
    df$dc_rel <- ifelse(df$dc_g > 0, df$dc_t - df$dc_g, -1000)
    df$dc_raw <- NULL

    unit_g <- df[!duplicated(df$dc_id), c("dc_id", "dc_g")]
    n_units <- nrow(unit_g)
    n_never <- sum(unit_g$dc_g == 0)
    cohorts <- table(unit_g$dc_g[unit_g$dc_g > 0])
    n_per <- length(unique(df$dc_t))
    balanced <- all(table(df$dc_id) == n_per)

    list(
        data = df,
        info = list(
            n_obs = nrow(df), n_units = n_units,
            n_treated = n_units - n_never, n_never = n_never,
            cohorts = cohorts,
            t_min = min(df$dc_t), t_max = max(df$dc_t), n_periods = n_per,
            balanced = balanced,
            n_missing = n_missing, n_always = n_always,
            covs = cov_map, cluster = !is.null(cluster), treat_type = treat_type,
            time_map = time_map, weighted = !is.null(weights), weight_label = weights,
            n_clusters = length(unique(df$dc_cl))))
}

# original period value for an internal one (see the recoding in dc_prepare)
dc_orig_time <- function(info, x) {
    tm <- info$time_map
    if (is.null(tm) || !isTRUE(tm$shifted)) x else tm$t0 + (x - 1) * tm$step
}

# HTML summary of a prepared panel (shared by all analyses)
dc_info_html <- function(i) {
    cohorts <- if (length(i$cohorts))
        paste(sprintf("%g (n=%d)", dc_orig_time(i, as.numeric(names(i$cohorts))), as.integer(i$cohorts)), collapse = ", ")
    else "none"
    msg <- c(
        sprintf("<b>%d</b> units, <b>%d</b> observations, periods %g to %g (%s panel).",
                i$n_units, i$n_obs, dc_orig_time(i, i$t_min), dc_orig_time(i, i$t_max),
                if (i$balanced) "balanced" else "unbalanced"),
        sprintf("Treated units: <b>%d</b>; never-treated units: <b>%d</b>.", i$n_treated, i$n_never),
        sprintf("Treatment cohorts (first treated period): %s.", cohorts))
    tm <- i$time_map
    if (!is.null(tm) && isTRUE(tm$shifted))
        msg <- c(msg, if (tm$step > 1)
            sprintf(paste0("Periods are %g apart: they are recoded to 1, 2, 3, ... for estimation, so event time e counts periods ",
                           "(e = 1 is %g time units after treatment). Periods and cohorts are shown with their original values."), tm$step, tm$step)
        else sprintf(paste0("Periods start at %g: they are recoded to start at 1 for estimation (0 is reserved for 'never treated'). ",
                            "Periods and cohorts are shown with their original values."), tm$t0))
    if (!is.null(tm) && isTRUE(tm$irregular))
        msg <- c(msg, paste0("Periods are not equally spaced: event time is measured in time units, and estimators that need the period ",
                             "just before treatment may fail or compare unequal gaps."))
    if (isTRUE(i$weighted))
        msg <- c(msg, sprintf("Sampling weights: %s (constant within units). Borusyak-Jaravel-Spiess is skipped with weights (didimputation has none).", i$weight_label))
    if (!is.null(i$n_clusters) && i$n_clusters < 30)
        msg <- c(msg, sprintf(paste0("<b>Only %d clusters.</b> Cluster-robust standard errors can be unreliable with few clusters; ",
                                     "Callaway-Sant'Anna's bootstrap and the t(G-1) intervals used for the regression estimators help only partly. ",
                                     "Treat small p-values with caution."), i$n_clusters))
    if (length(i$cohorts) && min(i$cohorts) < 5)
        msg <- c(msg, sprintf("The smallest treatment cohort has only %d unit(s); its effects are imprecise and a single unit can drive the aggregates.", min(i$cohorts)))
    if (!is.null(i$covs) && nrow(i$covs))
        msg <- c(msg, sprintf("Covariates: %s.", paste(i$covs$label, collapse = ", ")))
    if (!is.null(i$covs) && nrow(i$covs) && any(!i$covs$varying))
        msg <- c(msg, sprintf(paste0("Time-invariant covariates (%s) are absorbed by unit fixed effects; ",
                                     "only Callaway-Sant'Anna (and DDD) use them."),
                              paste(i$covs$label[!i$covs$varying], collapse = ", ")))
    if (!i$balanced) {
        msg <- c(msg, paste0("Unbalanced panel: Callaway-Sant'Anna, DDD and the pre-trend/HonestDiD tests use only units observed ",
                             "in every period; the fixed-effects estimators use all rows."))
        if (identical(i$treat_type, "indicator"))
            msg <- c(msg, paste0("Cohorts are inferred from the first observed treated period; if a unit's first treated row is ",
                                 "missing, supply the first treated period as the treatment variable instead."))
    }
    if (i$n_missing > 0)
        msg <- c(msg, sprintf("Dropped %d rows with missing or infinite values.", i$n_missing))
    if (i$n_always > 0)
        msg <- c(msg, sprintf("Dropped %d units already treated in the first period (no pre-period).", i$n_always))
    paste0("<p>", paste(msg, collapse = "<br>"), "</p>")
}

# keep a subset of rows and refresh the summary counts the estimators rely on
dc_subset <- function(d, keep) {
    df <- d$data[keep, , drop = FALSE]
    u <- df[!duplicated(df$dc_id), ]
    d$data <- df
    d$info$n_obs <- nrow(df); d$info$n_units <- nrow(u)
    d$info$n_never <- sum(u$dc_g == 0); d$info$n_treated <- sum(u$dc_g > 0)
    d$info$cohorts <- table(u$dc_g[u$dc_g > 0])
    d
}


# did's code for unbalanced panels (allow_unbalanced_panel = TRUE) is not safe to call here: on
# some inputs it crashes the R process. Panel estimators from the did family therefore run on
# the units observed in every period, and the number of units set aside is reported.
dc_balanced_subset <- function(d) {
    if (d$info$balanced)
        return(d)
    cnt <- table(d$data$dc_id)
    keep <- as.integer(names(cnt)[cnt == d$info$n_periods])
    if (!length(keep))
        stop("This estimator needs units observed in every period, and no unit has a complete record.")
    sub <- dc_subset(d, d$data$dc_id %in% keep)
    sub$info$balanced <- TRUE
    sub$info$n_unbalanced_dropped <- d$info$n_units - length(keep)
    sub
}
