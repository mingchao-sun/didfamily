# Core layer: "Data and setup guide" and "Method notes", shown at the top of each analysis when the
# matching options under "Notes and guide" are ticked (same layout as the rdfamily / scmfamily modules).
# Guides name the sample dataset, the exact field settings and the expected result; notes explain the
# method and every option on the analysis screen. Both are static, so they appear before any variable
# is chosen. Analyses that have run-time remarks append them under "This run".

dc_settings <- function(rows) {
    paste0("<table style='border-collapse:collapse;margin:4px 0 8px 0'>",
           "<tr><th style='text-align:left;padding:2px 14px 2px 0;border-bottom:1px solid #999'>Field</th>",
           "<th style='text-align:left;padding:2px 0;border-bottom:1px solid #999'>Set to</th></tr>",
           paste0("<tr><td style='padding:2px 14px 2px 0'>", names(rows), "</td><td style='padding:2px 0'>",
                  unlist(rows), "</td></tr>", collapse = ""), "</table>")
}

dc_layout <- function(cols, rows, extra = NULL) {
    th <- paste0("<th style='padding:2px 14px 2px 0;text-align:left;border-bottom:1px solid #999'>", cols, "</th>", collapse = "")
    tr <- paste0("<tr>", vapply(rows, function(r) paste0("<td style='padding:2px 14px 2px 0'>", r, "</td>", collapse = ""), ""),
                 "</tr>", collapse = "")
    paste0("<h4>Data layout</h4><p>One row per unit and period. Illustrative example:</p>",
           "<table style='border-collapse:collapse;margin:4px 0 8px 0'><tr>", th, "</tr>", tr, "</table>",
           "<ul>", extra, "</ul>")
}

# ---- guide ------------------------------------------------------------------------------------------

dc_panel_layout <- function(extra = NULL, cohort = FALSE)
    dc_layout(if (cohort) c("id", "year", "y", "g") else c("id", "year", "y", "D"),
              if (cohort) list(c("1", "2010", "4.2", "2012"), c("1", "2011", "4.5", "2012"), c("1", "2012", "5.9", "2012"),
                               c("2", "2010", "3.1", "0"), c("2", "2011", "3.3", "0"), c("2", "2012", "3.2", "0"))
              else list(c("1", "2010", "4.2", "0"), c("1", "2011", "4.5", "0"), c("1", "2012", "5.9", "1"),
                        c("2", "2010", "3.1", "0"), c("2", "2011", "3.3", "0"), c("2", "2012", "3.2", "0")),
              paste0("<li>Long panel: every unit appears once per period. Rows with a missing value in a variable used are dropped ",
                     "(the panel summary says how many). Units already treated in the first period are dropped because they have no pre-period.</li>",
                     "<li><b>Treatment variable</b>: a 0/1 indicator that turns on at the treatment date and never turns off, ",
                     "or the first treated period (0 or missing = never treated). Pick which under <i>Treatment variable is</i>.</li>",
                     "<li><b>Variable types in jamovi</b>: outcome Continuous; unit ID of type ID (or Nominal); time period Continuous with ",
                     "whole-number values (a year or a period number); treatment Continuous or Ordinal; covariates Continuous; ",
                     "cluster variable ID or Nominal; sampling weights Continuous and positive.</li>",
                     "<li>The cluster variable and the weights must be constant within a unit. Without a cluster variable, standard errors are clustered by unit.</li>",
                     extra))

dc_guide_staggered <- function(title, intro, expected, more = NULL) paste0(
    "<h3>", title, ": data and setup</h3><p>", intro, "</p>",
    "<h4>Sample data 1: simulated staggered adoption</h4>",
    dc_settings(list("Dataset" = "<b>Simulated staggered adoption</b>", "Outcome" = "y", "Unit ID" = "id",
                     "Time period" = "year", "Treatment variable" = "D", "Treatment variable is" = "A 0/1 indicator")),
    "<p>600 units, 12 periods, three cohorts (treated from period 5, 7 and 9) and 149 never-treated units. ", expected, "</p>",
    "<h4>Sample data 2: minimum wage and teen employment</h4>",
    dc_settings(list("Dataset" = "<b>Minimum wage and teen employment (mpdta, Callaway &amp; Sant'Anna)</b>", "Outcome" = "lemp", "Unit ID" = "countyreal",
                     "Time period" = "year", "Treatment variable" = "first.treat", "Treatment variable is" = "First treated period",
                     "Covariates" = "lpop (optional)")),
    "<p>500 counties, 2003 to 2007, cohorts 2004, 2006 and 2007. The 2004 to 2007 range gives a short event window: use Event window from -3 to 3.</p>",
    more,
    "<h4>Sample data 3: Castle Doctrine laws and homicide</h4>",
    dc_settings(list("Dataset" = "<b>Castle Doctrine laws and homicide (Cheng &amp; Hoekstra)</b>", "Outcome" = "l_homicide", "Unit ID" = "state",
                     "Time period" = "year", "Treatment variable" = "effyear", "Treatment variable is" = "First treated period",
                     "Sampling weights" = "pop2000 (optional)")),
    "<p>50 states, 2000 to 2010; cohorts 2005 to 2009 have only a few states each.</p>",
    dc_panel_layout())

dc_guide_raw <- function(key) {
    switch(key,

    didsetup = paste0(
        "<h3>Panel Structure Check: data and setup</h3>",
        "<p>Run this first. It reports whether the panel can be used by the estimators and what they will do with it.</p>",
        "<h4>Sample data</h4>",
        dc_settings(list("Dataset" = "<b>Simulated staggered adoption</b>", "Unit ID" = "id", "Time period" = "year",
                         "Treatment variable" = "D", "Outcome (optional)" = "y", "Treatment variable is" = "A 0/1 indicator")),
        "<p>Expected: 7200 usable rows, 600 units, 12 periods; balanced panel; treatment never turns off; 149 never-treated units (25%); ",
        "staggered adoption in 3 cohorts (5, 7, 9) and all checks OK. With <i>Castle Doctrine laws and homicide</i> (unit state, time year, treatment post) ",
        "the cohort table shows the few-units-per-cohort warning.</p>",
        dc_panel_layout("<li>The outcome is optional; with it the check also looks for missing and constant outcomes.</li>")),

    did22 = paste0(
        "<h3>Classic DiD (2x2): data and setup</h3>",
        "<p>One treated group, one comparison group, two periods. Gives regression DiD and the Sant'Anna-Zhao estimators in one table.</p>",
        "<h4>Sample data: simulated two-period DiD</h4>",
        dc_settings(list("Dataset" = "<b>Simulated two-period DiD</b>", "Outcome" = "y", "Treated group (0/1)" = "treated",
                         "Period (two values)" = "year", "Unit ID" = "id (optional)", "Covariates" = "x")),
        "<p>Expected with covariate x: regression DiD 1.82, doubly robust (improved) 1.48, inverse probability weighting 1.50, outcome regression 1.47. ",
        "The regression DiD is larger because it ignores the covariate imbalance between the groups.</p>",
        dc_layout(c("id", "year", "y", "treated", "x"),
                  list(c("1", "2010", "0.87", "0", "-0.9"), c("1", "2011", "1.31", "0", "-0.9"),
                       c("2", "2010", "1.02", "1", "0.4"), c("2", "2011", "3.20", "1", "0.4")),
                  paste0("<li>The treated-group variable is 0/1 and constant within a unit: 1 = the group that receives the treatment.</li>",
                         "<li>Exactly two period values are used (before and after). Without a unit ID, the data are treated as repeated cross-sections.</li>",
                         "<li><b>Variable types in jamovi</b>: outcome and covariates Continuous; treated group and period Nominal or Ordinal; unit ID ID or Nominal.</li>"))),

    didcompare = dc_guide_staggered("Compare DiD Estimators",
        "Runs up to six estimators on the same panel and puts their overall effects and event studies side by side.",
        "Expected overall ATT: TWFE 1.44, Callaway-Sant'Anna 1.64, Sun-Abraham 1.64, Borusyak-Jaravel-Spiess 1.61, Gardner 1.61, stacked 1.56. TWFE is lower because it uses already-treated units as controls while effects grow over time.",
        paste0("<p>Expected with covariate lpop, treatment as first treated period: Callaway-Sant'Anna about -0.042 (se 0.012). ",
               "With the castle data TWFE gives 0.08 and Callaway-Sant'Anna 0.11.</p>")),

    didtwfe = dc_guide_staggered("Two-Way Fixed Effects (TWFE)",
        "The conventional regression with unit and period fixed effects, with an event-study version.",
        "Expected: overall ATT 1.44 (se 0.034). For comparison the robust estimators give about 1.6.",
        "<p>Expected for the TWFE overall ATT without covariates: about -0.037.</p>"),

    didcs = dc_guide_staggered("Callaway & Sant'Anna",
        "Group-time effects ATT(g,t) aggregated to an overall effect and an event study.",
        "Expected (never-treated controls, no covariates): overall ATT 1.64 (se 0.051).",
        "<p>Expected with covariate lpop and the doubly robust method: overall ATT about -0.042 (se 0.012).</p>"),

    didsa = dc_guide_staggered("Sun & Abraham",
        "Interaction-weighted estimator with cohort-by-event-time effects; never-treated units are the controls.",
        "Expected: overall ATT 1.64 (se 0.038).",
        "<p>Sun and Abraham needs never-treated units; the panel check tells you whether there are any.</p>"),

    didbjs = dc_guide_staggered("Borusyak, Jaravel & Spiess (Imputation)",
        "Imputes each treated observation's untreated outcome from a model fitted on untreated observations only.",
        "Expected: overall ATT 1.61 (se 0.024).",
        "<p>This analysis has no sampling-weights field because the underlying package does not support weights.</p>"),

    didgardner = dc_guide_staggered("Gardner Two-Stage (did2s)",
        "Two-stage imputation estimator; point estimates equal the imputation estimator.",
        "Expected: overall ATT 1.61 (se 0.040); the same estimate as Borusyak-Jaravel-Spiess with a different standard error."),

    didetwfe = dc_guide_staggered("Extended TWFE (Wooldridge)",
        "Cohort-by-period treatment effects in a regression with cohort and period effects, with covariates, for linear, Poisson (counts) and logit (0/1) outcomes.",
        "Expected (linear model, not-yet-treated controls): overall ATT 1.61 (se 0.024), the same as the imputation estimators; with never-treated controls 1.64 (se 0.038) and pre-treatment effects -0.066 and -0.044 at e = -3 and -2.",
        paste0("<p>With covariate lpop and the linear model: overall ATT -0.051 (se 0.012), effects -0.033 at e = 0, -0.057 at e = 1 and -0.138 at e = 2.</p>",
               "<h4>Sample data 4: counts</h4>",
               dc_settings(list("Dataset" = "<b>Simulated staggered counts</b>", "Outcome" = "y", "Unit ID" = "id", "Time period" = "year", "Treatment variable" = "g",
                                "Treatment variable is" = "First treated period", "Covariates" = "x1", "Model for the outcome" = "Poisson")),
               "<p>600 units, 8 periods, cohorts 4, 5 and 6 and never-treated units, counts averaging 4.3. Expected (Poisson, not-yet-treated controls, Event window from -3 to 4): overall effect 2.60 additional events (se 0.144), by event time 1.76 (e = 0), 2.22 (e = 1), 2.78 (e = 2); by cohort 2.86 (cohort 4), 2.54 (5), 2.14 (6). ",
               "With never-treated controls: overall 2.51 (se 0.184) and pre-treatment effects -0.124 (e = -3) and 0.062 (e = -2), close to zero. With the linear model the overall effect is 2.72 (se 0.130).</p>")),
    didstacked = dc_guide_staggered("Stacked DiD",
        "One sub-experiment per cohort with clean controls, stacked and estimated with sub-experiment fixed effects.",
        "Expected: overall ATT 1.56 (se 0.029); 3 sub-experiments, 7795 stacked observations."),

    didiag = paste0(
        "<h3>Bacon, Pre-trend and HonestDiD: data and setup</h3>",
        "<p>Diagnoses how much to trust the TWFE estimate and how robust the event study is to violations of parallel trends.</p>",
        "<h4>Sample data</h4>",
        dc_settings(list("Dataset" = "<b>Simulated staggered adoption</b>", "Outcome" = "y", "Unit ID" = "id",
                         "Time period" = "year", "Treatment variable" = "D", "Goodman-Bacon decomposition" = "ticked",
                         "Joint pre-trend tests" = "ticked")),
        "<p>Expected: the Bacon decomposition reproduces the TWFE estimate 1.44 from nine 2x2 comparisons: treated vs untreated (weight 0.58, estimate 1.79), ",
        "earlier vs later treated (0.21, 0.59) and later vs earlier treated (0.21, 1.33). Joint pre-trend tests: TWFE rejects (chi-squared 75.9, p &lt; .001), ",
        "Sun-Abraham p = 0.36, Callaway-Sant'Anna p = 0.52. Tick <i>HonestDiD sensitivity analysis</i> only when you need it: it is slow.</p>",
        "<p><b>Pre-trend power</b> (Roth 2022): tick it in the section <i>Pre-trend power</i> (event study Sun &amp; Abraham, hypothesized trend 0.1). Expected: largest |t| among the pre-treatment effects 1.598 (the test passes), joint Wald p 0.362, ",
        "size 0.123 (with three pre-treatment effects a true-parallel-trends sample is rejected 12.3% of the time), power against a trend of 0.1 per period 1.000, and trends of 0.026 and 0.040 per period detected with 50% and 80% power. ",
        "The estimated effects are 0.611 at e = 0 and 1.157 at e = 1, so the trends that the test could miss are small compared with them.</p>",
        dc_panel_layout("<li>The Goodman-Bacon decomposition needs a balanced panel and a treatment that does not turn off.</li>")),

    dosedid = paste0(
        "<h3>Continuous Treatment DiD: data and setup</h3>",
        "<p>The treatment has an intensity (a dose). Estimates the dose-response ATT(d) and the marginal effect ACRT(d).</p>",
        "<h4>Sample data: simulated continuous treatment</h4>",
        dc_settings(list("Dataset" = "<b>Simulated continuous treatment</b>", "Outcome" = "y", "Unit ID" = "id",
                         "Time period" = "year", "First treated period" = "g", "Dose" = "dose")),
        "<p>Expected with the default linear spline: average effect (ATT) 1.44 (se 0.043), average marginal effect of dose (ACRT) 3.10 (se 0.037).</p>",
        dc_layout(c("id", "year", "y", "g", "dose"),
                  list(c("1", "1", "-0.28", "3", "0.21"), c("1", "3", "0.90", "3", "0.21"), c("2", "1", "0.10", "0", "0"), c("2", "3", "0.05", "0", "0")),
                  paste0("<li><b>First treated period</b>: 0 = never treated. <b>Dose</b>: constant within a unit, positive for treated units and 0 for the untreated.</li>",
                         "<li>The estimator needs units observed in every period; units with gaps are set aside and counted in the notes.</li>",
                         "<li><b>Variable types in jamovi</b>: outcome and dose Continuous; unit ID ID; time period Continuous (whole numbers).</li>"))),

    dddcompare = paste0(
        "<h3>Triple Differences (DDD): data and setup</h3>",
        "<p>Compares eligible with ineligible units, in treated versus untreated groups, to remove a shock that hits the whole treated group.</p>",
        "<h4>Sample data: simulated triple differences</h4>",
        dc_settings(list("Dataset" = "<b>Simulated triple differences</b>", "Outcome" = "y", "Unit ID" = "id", "Time period" = "year",
                         "First treated period" = "g", "Eligibility (0/1)" = "elig")),
        "<p>Expected: triple difference 1.51 (se 0.059); the DiD for eligible units alone is 2.32 and for ineligible units alone 0.80, so a group-wide shock of about 0.8 is removed.</p>",
        dc_layout(c("id", "year", "y", "g", "elig"),
                  list(c("1", "1", "-0.74", "4", "0"), c("1", "4", "1.20", "4", "0"), c("2", "1", "0.31", "4", "1"), c("2", "4", "3.05", "4", "1")),
                  paste0("<li><b>First treated period</b> belongs to the group (0 = never-treated group); <b>Eligibility</b> is 0/1 within every group, constant within a unit.</li>",
                         "<li>Both eligible and ineligible units are needed in treated and in untreated groups.</li>"))),

    dcdh = paste0(
        "<h3>Treatment Switching On and Off: data and setup</h3>",
        "<p>Dynamic effects of a binary treatment that can turn on and off again, with placebo tests (de Chaisemartin and D'Haultfoeuille, 2024). The other analyses need a treatment that never turns off.</p>",
        "<h4>Sample data: simulated switching treatment</h4>",
        dc_settings(list("Dataset" = "<b>Simulated treatment that switches on and off</b>", "Outcome" = "y", "Group ID" = "id", "Time period" = "year",
                         "Treatment (0/1)" = "D", "Number of effects" = "3", "Number of placebos" = "2")),
        "<p>300 groups, 2001 to 2010. 218 groups switch the treatment on between 2004 and 2007 and about a third of them switch it off again two periods later; 82 never switch. ",
        "The effect of the treatment is 2. Expected: effect 1 = 1.995 (s.e. 0.108), effect 2 = 1.891 (0.119), effect 3 = 1.168 (0.134), which falls because part of the switchers have left the treatment by then; ",
        "placebos -0.096 (0.112) and 0.026 (0.113), joint placebo test p = 0.57.</p>",
        dc_layout(c("id", "year", "y", "D"),
                  list(c("1", "2004", "6.1", "0"), c("1", "2005", "8.4", "1"), c("1", "2006", "8.0", "1"), c("1", "2007", "6.3", "0")),
                  paste0("<li>Long panel: one row per group and period. The treatment is 0/1 and <b>may switch on and off</b>; a nominal treatment with two values is read with the second value as treated.</li>",
                         "<li>Periods are whole numbers. Missing periods are allowed (groups are compared when observed); rows with missing values are dropped.</li>",
                         "<li>Needed: groups that change the treatment, and comparison groups with the same starting treatment that have not changed it yet.</li>"))),

    "")
}

dc_guide_html <- function(key) {
    x <- dc_guide_raw(key)
    sub("</p>", paste0("</p><p><i>Open the sample data from File &gt; Open &gt; Data Library &gt; Difference-in-Differences Family; ",
                       "the library lists each data set by the name shown below.</i></p>"), x, fixed = TRUE)
}

# ---- method notes -----------------------------------------------------------------------------------

dc_n_vars <- function(opt = NULL, weights = TRUE, cluster = TRUE, covs = TRUE) paste0(
    "<h2>Data</h2>",
    "<p>The data must be a <b>long panel</b>: one row per unit and period. Rows with missing values in the variables used are dropped; ",
    "units treated in the first period are dropped because they have no pre-period.</p><ul>",
    "<li><b>Outcome</b>: numeric variable to analyse.</li>",
    "<li><b>Unit ID</b>: identifier of the unit (for example county or state).</li>",
    "<li><b>Time period</b>: whole-number period (year or period index).</li>",
    "<li><b>Treatment variable</b>: a 0/1 indicator that switches on and stays on, or the first treated period (0 or missing = never treated). ",
    "<b>Treatment variable is</b> tells the module which one you supplied.</li>",
    if (covs) "<li><b>Covariates</b> (optional): variables that make the comparison more credible. How they enter depends on the estimator (see below).</li>",
    if (cluster) "<li><b>Cluster variable</b> (optional): standard errors are clustered on it; it must be constant within a unit. The default is the unit.</li>",
    if (weights) "<li><b>Sampling weights</b> (optional): positive, constant within a unit.</li>",
    opt, "</ul>")

dc_n_event <- function() paste0(
    "<h2>Event window and confidence level</h2><ul>",
    "<li><b>Event window, from / to</b>: the leads and lags shown in the dynamic table and the plot. Event time e = 0 is the first treated period; ",
    "the reference period is e = -1 (its effect is fixed at zero). The window must contain e = -1 and e = 0.</li>",
    "<li><b>Confidence level (%)</b>: for the intervals and p-values.</li></ul>")

dc_n_plot <- function(dodge = FALSE) paste0(
    "<h2>Plot style</h2><ul><li><b>Plot colors</b>: color or gray scale (for print).</li>",
    "<li><b>Show confidence intervals as</b>: error bars or a shaded ribbon.</li>",
    if (dodge) "<li><b>Spacing between estimators</b>: horizontal offset of the error bars so they do not overlap.</li>",
    "</ul>")

dc_n_code <- function() paste0(
    "<li><b>Show reproducible R code</b>: prints a script that repeats the estimation outside jamovi with the same settings.</li>")

dc_n_reading <- function(...) paste0("<h2>Reading the results</h2><ul>", paste0("<li>", c(...), "</li>", collapse = ""), "</ul>")

dc_n_inference <- function(t = FALSE) paste0(
    "<li><b>Inference</b>: standard errors are clustered at the unit (or the cluster variable); ",
    if (t) "confidence intervals and p-values use t(G-1) critical values with G clusters, as fixest and Stata's reghdfe do. "
    else "confidence intervals and p-values use the normal approximation, as in the underlying package. ",
    "With few clusters (the panel summary warns you) they can be too optimistic.</li>")

dc_method_html <- function(key, o = NULL) {
    H <- function(title, intro) paste0("<h1>", title, "</h1><p>", intro, "</p>")
    single <- function(k, title, intro, est, extra_opts = NULL, limits, wts = TRUE, covs_txt = NULL, tinf = FALSE) paste0(
        H(title, intro),
        dc_n_vars(weights = wts, covs = TRUE),
        "<h2>Estimation</h2><ul>", est, "</ul>",
        extra_opts,
        dc_n_event(),
        "<h2>Additional output</h2><ul>", dc_n_code(), "</ul>",
        dc_n_plot(),
        dc_n_reading(
            "<b>Overall ATT</b>: the average effect on treated units after treatment, with its standard error, confidence interval and p-value.",
            "<b>Dynamic effects</b>: one effect per event time. Leads (e &lt; 0) should be close to zero if parallel trends is plausible; lags (e &ge; 0) show how the effect evolves.",
            "<b>Event study plot</b>: the dynamic effects with their intervals; the reference period -1 is fixed at zero.",
            dc_n_inference(tinf)),
        "<h2>Notes and limitations</h2><ul>", limits,
        "<li>To compare this estimator with the others on the same panel, use <i>Compare DiD Estimators</i>.</li></ul>")
    switch(key,

    didsetup = paste0(
        H("Panel Structure Check", "Describes the panel before any estimation: size, balance, duplicates, treatment timing and cohorts."),
        dc_n_vars(weights = FALSE, cluster = FALSE, covs = FALSE, opt = "<li><b>Outcome</b> is optional here: with it, the check also looks for missing and constant outcomes.</li>"),
        "<h2>Plot style</h2><ul><li><b>Plot colors</b>: color or gray scale.</li></ul>",
        dc_n_reading(
            "<b>Panel checks</b>: <i>Problem</i> means the estimators will refuse the data or give wrong answers; <i>Check</i> means something will be dropped or inferred; <i>Note</i> is informational.",
            "<b>Treatment cohorts</b>: units per cohort and periods before and after treatment. Few units per cohort make cohort-specific estimates noisy.",
            "<b>Treatment status plot</b>: each row is a unit and each column a period. With staggered adoption the treated cells form a staircase; a single block means one treatment date. With many units an evenly spaced sample is drawn."),
        "<h2>Notes and limitations</h2><ul>",
        "<li>No never-treated units, or very short pre-periods, limit which estimators are reliable (Sun-Abraham needs never-treated units).</li>",
        "<li>Run this first, then pick an estimator from the Estimators menu or use <i>Compare DiD Estimators</i>.</li></ul>"),

    did22 = paste0(
        H("Classic DiD (2x2)", "Compares the change in the outcome of the treated group between two periods with the change in the comparison group. The estimate is the average effect on the treated (ATT); it needs parallel trends (without the treatment the groups would have changed alike) and no anticipation."),
        "<h2>Data</h2><ul>",
        "<li><b>Outcome</b>: numeric variable.</li><li><b>Treated group (0/1)</b>: 1 = the group that receives the treatment.</li>",
        "<li><b>Period (two values)</b>: before and after.</li>",
        "<li><b>Unit ID</b> (optional): with it the data are treated as a panel; without it as repeated cross-sections.</li>",
        "<li><b>Covariates</b> (optional): used by the doubly robust, inverse probability weighting and outcome-regression estimators; regression DiD enters them linearly.</li></ul>",
        "<h2>Estimators</h2><ul>",
        "<li><b>Regression DiD (OLS)</b>: the outcome regressed on group, period and their interaction; the interaction is the ATT and equals the difference in changes of the group means when there are no covariates.</li>",
        "<li><b>Doubly robust (improved)</b>: Sant'Anna and Zhao's locally efficient estimator; it is consistent if the propensity-score model or the outcome model is right.</li>",
        "<li><b>Doubly robust (traditional)</b>: the earlier doubly robust form, without the efficiency improvement.</li>",
        "<li><b>Inverse probability weighting</b>: reweights the comparison group to resemble the treated group on the covariates, using a propensity score.</li>",
        "<li><b>Outcome regression</b>: predicts the treated group's change from the comparison group's covariate relationship.</li></ul>",
        "<h2>Confidence level and plot</h2><ul><li><b>Confidence level (%)</b>: for the intervals and p-values.</li><li><b>Plot colors</b>: color or gray scale.</li></ul>",
        dc_n_reading(
            "<b>Mean outcome by group and period</b>: the four cell means behind the 2x2 comparison.",
            "<b>DiD estimates</b>: the ATT from each selected estimator with its standard error, interval and p-value. Without covariates all estimators coincide.",
            "<b>Group means plot</b>: the two groups before and after; parallel lines would indicate no effect."),
        "<h2>Notes and limitations</h2><ul>",
        "<li>Two periods cannot test parallel trends; use the staggered-adoption analyses with several pre-periods for that.</li>",
        "<li>If the estimates differ a lot between regression DiD and the covariate-adjusted ones, the groups differ in covariates that matter.</li></ul>"),

    didcompare = paste0(
        H("Compare DiD Estimators", "Runs the heterogeneity-robust estimators and TWFE on the same panel so that overall effects and event studies can be compared directly."),
        dc_n_vars(opt = "<li>Sampling weights are not used by Borusyak-Jaravel-Spiess (the underlying package has none); that row is skipped when weights are given.</li>"),
        "<h2>Estimators</h2><p>Tick the estimators to include. What each does:</p><ul>",
        paste0("<li>", vapply(c("twfe", "cs", "sa", "bjs", "gardner", "stacked"), function(k) dc_single_notes(k, o), ""), "</li>", collapse = ""),
        "</ul>",
        "<h2>Callaway &amp; Sant'Anna settings</h2><ul>",
        "<li><b>Control group</b>: never-treated units, or not-yet-treated units.</li>",
        "<li><b>Estimation method (with covariates)</b>: doubly robust, inverse probability weighting or outcome regression.</li>",
        "<li><b>Bootstrap standard errors</b>, <b>Bootstrap draws</b>, <b>Random seed</b>: multiplier bootstrap in place of analytical standard errors; the seed makes it reproducible.</li></ul>",
        dc_n_event(),
        "<h2>Additional output</h2><ul>",
        "<li><b>Forest plot of overall ATT</b>: the overall effects with intervals, one row per estimator.</li>",
        "<li><b>Group-time effects (heat map)</b>, <b>Effects by cohort</b>, <b>Effects by calendar time</b>: Callaway &amp; Sant'Anna aggregations.</li>",
        dc_n_code(), "</ul>",
        dc_n_plot(dodge = TRUE),
        dc_n_reading(
            "<b>Overall ATT by estimator</b>: the average post-treatment effect from each estimator. Large gaps between TWFE and the robust estimators point to heterogeneous or dynamic effects.",
            "<b>Dynamic effects</b> and the <b>event study comparison</b>: all estimators on one scale; leads should be near zero.",
            "Estimators that cannot run (for example Sun-Abraham without never-treated units) show a note instead of numbers."),
        "<h2>Notes and limitations</h2><ul>",
        "<li>The estimators have different control groups and weighting, so overall ATTs need not agree even when all are valid.</li>",
        "<li>Each estimator also has its own analysis in the Estimators menu, with identical numbers.</li></ul>"),

    didtwfe = single("twfe", "Two-Way Fixed Effects (TWFE)", dc_single_notes("twfe", o),
        "<li>The static estimate is the coefficient on the treatment indicator in a regression with unit and period fixed effects. The event study replaces it by leads and lags with e = -1 omitted. Covariates enter linearly.</li>",
        limits = "<li>With staggered adoption and effects that differ across cohorts or grow over time, TWFE can put negative weights on some comparisons and differ from the robust estimators even in sign. <i>Bacon, Pre-trend and HonestDiD</i> shows the weights.</li>", tinf = TRUE),

    didcs = single("cs", "Callaway & Sant'Anna", dc_single_notes("cs", o),
        paste0("<li><b>Control group</b>: <i>Never treated</i> uses only units never treated; <i>Not yet treated</i> also uses units treated later, and needs fewer never-treated units.</li>",
               "<li><b>Estimation method (with covariates)</b>: doubly robust (default), inverse probability weighting or outcome regression.</li>",
               "<li><b>Multiplier bootstrap</b>, <b>Bootstrap draws</b>, <b>Random seed</b>: bootstrap standard errors in place of analytical ones. It is used automatically when clustering above the unit level.</li>"),
        extra_opts = paste0("<h2>Additional output</h2><ul>",
            "<li><b>Group-time effects (heat map)</b>: ATT(g,t) for every cohort g and period t.</li>",
            "<li><b>Effects by cohort</b>: the overall effect for each treatment cohort.</li>",
            "<li><b>Effects by calendar time</b>: the average effect among treated units in each period.</li></ul>"),
        limits = "<li>Universal base period (g-1). On unbalanced panels only units observed in every period are used. Cohorts with very few units give noisy ATT(g,t).</li>"),

    didsa = single("sa", "Sun & Abraham", dc_single_notes("sa", o),
        "<li>The regression interacts cohort with event-time dummies, which removes the contamination of TWFE. Never-treated units are the controls, so the analysis needs some. Covariates that vary within units enter linearly.</li>",
        limits = "<li>Without never-treated units the estimator cannot run (the panel summary shows how many there are).</li>", tinf = TRUE),

    didbjs = single("bjs", "Borusyak, Jaravel & Spiess (Imputation)", dc_single_notes("bjs", o),
        "<li>Unit and period effects (and covariate coefficients) are fitted on untreated observations only; each treated observation's counterfactual is imputed and the effect is the average treated-minus-imputed gap.</li>",
        limits = "<li>Sampling weights are not available. Pre-treatment coefficients come from a separate pre-trend regression on untreated observations, so they are a real test of the assumption.</li>", wts = FALSE),

    didgardner = single("gardner", "Gardner Two-Stage (did2s)", dc_single_notes("gardner", o),
        "<li>Stage 1 fits unit and period effects on untreated observations; stage 2 regresses the adjusted outcome on the treatment. Point estimates equal the imputation estimator; the variance estimator accounts for the first stage.</li>",
        limits = "<li>Pre-treatment coefficients are close to zero by construction (the first stage also uses the pre-periods), so they are not a pre-trend test.</li>"),

    didetwfe = single("etwfe", "Extended TWFE (Wooldridge)", dc_single_notes("etwfe", o),
        paste0("<li><b>Control group</b>: <i>Not yet treated</i> (default) uses units not yet treated as well as never-treated units; <i>Never treated</i> uses only never-treated units and also gives pre-treatment effects (leads) for the pre-trend check.</li>",
               "<li><b>Model for the outcome</b>: linear for any outcome; Poisson for counts or other non-negative outcomes (effects are in the units of the outcome, not log points); logit for 0/1 outcomes (effects are in probability). The nonlinear models have no fixed effects, they use cohort and period dummies.</li>",
               "<li>Covariates are demeaned within cohort and interacted with the cell effects, so the effect can vary with them; the overall effect averages these effects over the treated observations.</li>"),
        extra_opts = paste0("<h2>Additional output</h2><ul>",
            "<li><b>Effects by cohort</b> and <b>Effects by calendar time</b>: the same effects averaged over each cohort or each period.</li></ul>"),
        limits = "<li>No sampling weights. The overall effect weights observations equally, so it can differ from estimators that weight cohorts differently. Without covariates and with a linear model the point estimates equal the imputation estimators. The Poisson and logit standard errors are clustered and delta-method based.</li>", wts = FALSE, tinf = TRUE),
    didstacked = single("stacked", "Stacked DiD", dc_single_notes("stacked", o),
        "<li>One sub-experiment per cohort, with that cohort and its clean controls, stacked with unit-by-experiment and period-by-experiment fixed effects. The overall effect weights cohorts by size and by how many controls each uses.</li>",
        limits = "<li>The overall ATT need not equal those of the other estimators because of the weighting.</li>", tinf = TRUE),

    didiag = paste0(
        H("Bacon, Pre-trend and HonestDiD", "Diagnostics for staggered adoption: what the TWFE estimate is made of, whether pre-treatment effects are jointly zero, and how robust the effect is to violations of parallel trends."),
        dc_n_vars(weights = FALSE, cluster = FALSE, covs = FALSE),
        "<h2>Diagnostics</h2><ul>",
        "<li><b>Goodman-Bacon decomposition</b>: the TWFE estimate is a weighted average of all 2x2 comparisons. Comparisons of later-treated with earlier-treated units use already-treated units as controls and can bias the estimate. Needs a balanced panel.</li>",
        "<li><b>TWFE weights on treated observations</b>: how much each treated observation counts; negative weights mean TWFE can differ from robust estimators even in sign.</li>",
        "<li><b>Joint pre-trend tests</b>: Wald tests that all pre-treatment effects are zero, for TWFE, Sun-Abraham and Callaway-Sant'Anna. Failing to reject is not proof of parallel trends because the test can have low power.</li></ul>",
        "<h2>Sensitivity analysis</h2><ul>",
        "<li><b>HonestDiD sensitivity analysis</b> (slow, off by default): confidence intervals that stay valid when post-treatment violations of parallel trends are up to Mbar times the largest pre-treatment violation (relative magnitudes).</li>",
        "<li><b>Event study to analyse</b>: Callaway-Sant'Anna or Sun-Abraham.</li>",
        "<h2>Pre-trend power</h2><ul><li><b>Power of the pre-trend test</b> (Roth 2022): the usual pre-test passes if no pre-treatment effect is individually significant. The table gives the probability that it rejects when trends are parallel (size, above the nominal level when there are several pre-treatment effects), the probability that it rejects against a linear trend you specify (<i>Hypothesized trend</i>, slope per period through the reference period e = -1), and the slopes it detects with 50 and 80 percent power. A second table shows the bias such trends would put in each post-treatment effect, to compare with the estimate; the plot draws the trends over the event study.</li></ul>",
        "<li><b>Effect of interest</b>: the first post-treatment effect, or the average of the post-treatment effects.</li>",
        "<li><b>Largest relative magnitude (Mbar)</b>: 0 reproduces the original interval; larger values allow bigger violations. The plot shows how the interval grows with Mbar.</li></ul>",
        "<h2>Event window and settings</h2><ul>",
        "<li><b>Callaway &amp; Sant'Anna control group</b>: never treated or not yet treated.</li>",
        "<li><b>Event window</b> and <b>Confidence level (%)</b>: as in the estimators.</li></ul>",
        "<h2>Plot style</h2><ul><li><b>Plot colors</b>: color or gray scale.</li></ul>",
        dc_n_reading(
            "A TWFE estimate that differs from the robust estimators, together with a large weight on later-vs-earlier comparisons or negative weights, indicates contamination.",
            "A joint test that rejects for TWFE but not for the robust estimators suggests that the TWFE leads reflect contamination, not a pre-trend.",
            "The HonestDiD breakdown value is the smallest Mbar at which the interval includes zero; a larger value means a more robust result."),
        "<h2>Notes and limitations</h2><ul><li>HonestDiD is run on the event study, so its result depends on the window and the chosen estimator.</li></ul>"),

    dosedid = paste0(
        H("Continuous Treatment DiD", "Callaway, Goodman-Bacon and Sant'Anna (2024). Estimates the dose-response ATT(d), the effect of dose d relative to no treatment, and the marginal effect ACRT(d) of one more unit of dose, with their averages."),
        "<h2>Data</h2><ul>",
        "<li><b>Outcome</b>, <b>Unit ID</b>, <b>Time period</b>: as in the other analyses; the panel must be observed in every period.</li>",
        "<li><b>First treated period</b>: 0 = never treated.</li>",
        "<li><b>Dose</b>: treatment intensity, constant within a unit; positive for treated units and 0 for the untreated.</li></ul>",
        "<h2>Estimation</h2><ul>",
        "<li><b>Control group</b>: never-treated units or not-yet-treated units.</li>",
        "<li><b>Spline degree</b> and <b>Number of interior knots</b>: the dose-response is a B-spline of this degree. Degree 1 with no knots is a straight line; add knots or raise the degree for curvature.</li>",
        "<li><b>Confidence level (%)</b> and <b>Plot colors</b>.</li></ul>",
        dc_n_reading(
            "<b>Average effects over all treated units</b>: the overall ATT and the average marginal effect ACRT.",
            "<b>Dose-response at selected doses</b>: ATT(d) and ACRT(d) at doses across the observed range.",
            "<b>Dose-response curves</b>: the same as plots; pointwise intervals, not simultaneous bands."),
        "<h2>Notes and limitations</h2><ul>",
        "<li>ATT(d) needs parallel trends for untreated outcomes. Comparing ATT(d) across doses, and reading ACRT(d) as a causal marginal effect, needs the stronger assumption that there is no selection into dose on the size of the effect.</li>",
        "<li>Covariates are not used here. If the curve looks off near the smallest doses, raise the degree or add knots and check that the shape is stable.</li></ul>"),

    dddcompare = paste0(
        H("Triple Differences (DDD)", "Removes a shock that hits the whole treated group, eligible or not: DDD is the DiD for eligible units minus the DiD for ineligible units."),
        dc_n_vars(weights = FALSE, cluster = FALSE, opt = "<li><b>First treated period</b> belongs to the group (0 = never-treated group); <b>Eligibility (0/1)</b>: 1 = eligible for the policy, constant within a unit.</li>"),
        "<h2>Estimation</h2><ul>",
        "<li><b>Control group</b>: never-treated groups or not-yet-treated groups.</li>",
        "<li><b>Estimation method</b>: doubly robust, inverse probability weighting or outcome regression (they differ only when covariates are used).</li></ul>",
        dc_n_event(), dc_n_plot(dodge = TRUE),
        dc_n_reading(
            "<b>Overall ATT</b>: the triple-difference estimate next to the DiD for eligible units and the DiD for ineligible units. The ineligible DiD shows the size of the group-wide shock.",
            "<b>Dynamic effects</b> and the <b>event study comparison</b>: the three series by event time; the DDD leads should be near zero."),
        "<h2>Notes and limitations</h2><ul>",
        "<li>DDD needs parallel trends in the <i>difference</i> between eligible and ineligible units, not within each. It can hold even when the single DiDs fail.</li>",
        "<li>Both eligible and ineligible units are needed in treated and untreated groups.</li></ul>"),

    dcdh = paste0(
        H("Treatment Switching On and Off", "de Chaisemartin and D'Haultfoeuille (2024). Estimates the effect of a binary treatment when groups can switch it on and off, which two-way fixed effects and the staggered-adoption estimators do not allow."),
        "<h2>Data</h2><ul>",
        "<li><b>Outcome</b>, <b>Group ID</b>, <b>Time period</b>: a long panel with one row per group and period; periods are whole numbers.</li>",
        "<li><b>Treatment (0/1)</b>: may turn on and off. A group's <i>starting treatment</i> is its value in its first period; its <i>first change</i> is the first period with a different value.</li></ul>",
        "<h2>Estimation</h2><ul>",
        "<li><b>Effect l</b>: for every group that changes its treatment, the change in the outcome from the period before the first change to l-1 periods after it, minus the same change for groups with the same starting treatment that had not changed it by then; averaged over those groups. ",
        "Effect 1 is the change in the period of the first switch. Groups that switch off are counted with the sign reversed, so every effect is the effect of the treatment relative to no change.</li>",
        "<li><b>Placebo l</b>: the same comparison l periods before the period preceding the first change. It should be near zero if groups that switch would have followed the same trend as groups that do not.</li>",
        "<li><b>Number of effects</b> and <b>Number of placebos</b>: how many are estimated; later ones use fewer groups.</li>",
        "<li><b>Confidence level (%)</b>: for the intervals; the standard errors are the influence-function estimates of the paper, which allow for the estimated comparison means.</li></ul>",
        "<h2>Plot style</h2><ul><li><b>Event study plot</b> shows the placebos (negative distances) and the effects; <b>Plot colors</b>: color or gray scale.</li></ul>",
        dc_n_reading(
            "<b>Effects</b>: the average effect of having switched, k periods after the switch. Because some groups switch back, later effects mix groups that are still treated and groups that left the treatment, so they need not equal the effect of a lasting treatment.",
            "<b>Placebos</b> and the joint test: a rejection signals that groups which switch were already on a different path.",
            "<b>Switchers</b>: the number of groups used by that estimate."),
        "<h2>Notes and limitations</h2><ul>",
        "<li>Binary treatment, no covariates, no weights in this version; use the estimators in this menu for staggered adoption with an absorbing treatment.</li>",
        "<li>Results equal the R package DIDmultiplegtDYN and the Stata command did_multiplegt_dyn on the same data (balanced and unbalanced panels).</li></ul>"),

    "")
}

# static method notes followed by remarks that depend on the run
dc_run_notes <- function(key, o, run_html)
    paste0(dc_method_html(key, o), "<h2>This run</h2>", run_html)

# Tables, plots and code blocks stay hidden until the required variables are assigned, so that a fresh
# analysis shows no empty tables (as in rdfamily / scmfamily). Once all are assigned the visibility
# rules of the .r.yaml (showCode, aggGroup, ...) apply unchanged.
dc_hide_empty <- function(self, req) {
    o <- self$options
    if (all(vapply(req, function(k) length(o[[k]]) > 0, TRUE))) return(invisible())
    for (nm in c("means", "summary", "dynamic", "plot", "forest", "groupTab", "calTab", "gtPlot", "code", "curve",
                 "checks", "cohorts", "bacon", "baconPlot", "weights", "pretest", "honest", "honestPlot")) {
        el <- tryCatch(self$results[[nm]], error = function(e) NULL)
        if (!is.null(el)) el$setVisible(FALSE)
    }
    invisible()
}
