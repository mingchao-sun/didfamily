# Core layer: DiD with a continuous treatment (Callaway, Goodman-Bacon & Sant'Anna 2024) via contdid.
# Needs a first-treated-period variable (0 = never treated) and a dose that is constant within
# a unit (0 for untreated units).

# Returns list(ok, message, att, att_se, acrt, acrt_se, curve) where `curve` has one row per dose
# on contdid's grid: dose, att, att_se, acrt, acrt_se (pointwise standard errors).
dc_dose <- function(d, control = c("nevertreated", "notyettreated"), degree = 1, knots = 0) {
    control <- match.arg(control)
    # contdid draws a multiplier bootstrap for the overall standard errors even with bstrap = FALSE:
    # fix the seed so that the result is reproducible, and put the caller's random number state back afterwards
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    set.seed(1234)
    tryCatch({
        d <- dc_balanced_subset(d)
        if (control == "nevertreated" && d$info$n_never == 0)
            stop("No never-treated units: choose 'Not yet treated' as control group.")
        fit <- contdid::cont_did(
            yname = "dc_y", dname = "dc_dose", gname = "dc_g", tname = "dc_t", idname = "dc_id",
            data = d$data, target_parameter = "level", aggregation = "dose",
            treatment_type = "continuous", dose_est_method = "parametric",
            degree = as.integer(degree), num_knots = as.integer(knots),
            control_group = control, bstrap = FALSE, cband = FALSE, print_details = FALSE)
        curve <- data.frame(dose = as.numeric(fit$dose), att = as.numeric(fit$att.d),
                            att_se = as.numeric(fit$att.d_se), acrt = as.numeric(fit$acrt.d),
                            acrt_se = as.numeric(fit$acrt.d_se))
        list(ok = TRUE, message = "", att = fit$overall_att, att_se = fit$overall_att_se,
             acrt = fit$overall_acrt, acrt_se = fit$overall_acrt_se, curve = curve,
             dropped = d$info$n_unbalanced_dropped)
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
}
