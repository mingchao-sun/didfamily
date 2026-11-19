# Core layer: triple differences (DDD). Needs a first-treated-period variable defined for
# the unit's group and a 0/1 eligibility (partition) variable within the group.
# Only eligible units in treated groups are affected by the policy; ineligible units share
# group-level shocks, which the third difference removes.

dc_ddd <- function(d, min_e, max_e, control = c("nevertreated", "notyettreated"),
                   est_method = c("dr", "ipw", "reg")) {
    method <- "Triple difference (DDD)"
    control <- match.arg(control)
    est_method <- match.arg(est_method)
    tryCatch({
        d <- dc_balanced_subset(d)
        if (control == "nevertreated" && d$info$n_never == 0)
            stop("No never-treated groups: choose 'Not yet treated' as control group.")
        cx <- dc_cov_expr(d, fe = FALSE)
        xf <- if (is.null(cx)) ~1 else eval(bquote(~ .(cx)))
        fit <- triplediff::ddd(yname = "dc_y", tname = "dc_t", idname = "dc_id", gname = "dc_g",
                               pname = "dc_p", xformla = xf, data = d$data,
                               control_group = control, base_period = "universal",
                               est_method = est_method)
        simple <- triplediff::agg_ddd(fit, type = "simple")$aggte_ddd
        es <- triplediff::agg_ddd(fit, type = "eventstudy", min_e = min_e, max_e = max_e)$aggte_ddd
        dc_out(method, att = simple$overall.att, se = simple$overall.se,
               dynamic = dc_dyn(es$egt, es$att.egt, es$se.egt, min_e, max_e),
               control = if (control == "nevertreated") "Never-treated groups" else "Not-yet-treated groups")
    }, error = function(e) dc_fail(method, e))
}

# Ordinary Callaway-Sant'Anna DiD within one eligibility group, for comparison with the DDD.
dc_did_part <- function(d, which = c("eligible", "ineligible"), min_e, max_e, control, est_method) {
    which <- match.arg(which)
    method <- if (which == "eligible") "DiD, eligible only" else "DiD, ineligible only"
    sub <- dc_subset(d, d$data$dc_p == if (which == "eligible") 1 else 0)
    r <- dc_cs(sub, min_e, max_e, control, est_method)
    r$method <- method
    r
}
