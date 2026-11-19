# didfamily — compare difference-in-differences estimators in jamovi

A jamovi module for running several modern DiD estimators on the same panel and seeing where they
agree and disagree, with diagnostics, sensitivity analysis and reproducible R code.
Version 0.5.0 (development). Chinese summary: [README.zh.md](README.zh.md).

## Analyses (jamovi menu: DiD)

The menu is in workflow order, with sub-headings: panel check and 2x2 first, then single estimators, then comparison, then diagnostics and extensions.

| Analysis | What it does |
|---|---|
| **Panel Structure Check** | Checks the panel before estimating (balance, duplicates, gaps, missing values, absorbing treatment, never-treated share, cohort sizes), cohort table, treatment-status plot. Problems are listed, not thrown. |
| **Classic DiD (2x2)** | Two groups, two periods, with or without a unit ID: group means, regression DiD, doubly robust / IPW / outcome-regression estimators (`DRDID`), plot with the implied counterfactual. |
| *Estimators* | One analysis per estimator, each with its own options, overall ATT, dynamic effects table, event-study plot, optional R code: **Two-Way Fixed Effects (TWFE)**, **Callaway & Sant'Anna** (control group, doubly robust / IPW / regression, multiplier bootstrap, group-time heat map, cohort and calendar-time effects), **Sun & Abraham**, **Borusyak, Jaravel & Spiess (Imputation)**, **Gardner Two-Stage (did2s)**, **Extended TWFE (Wooldridge)** (linear, Poisson and logit models, covariates, not-yet or never-treated controls), **Stacked DiD**. They call the same code as the comparison, so the numbers agree exactly. |
| **Compare DiD Estimators** | Up to six estimators on one staggered panel (table below), overall ATT table, dynamic effects, overlaid event-study plot, forest plot, Callaway-Sant'Anna group-time heat map and cohort / calendar-time effects, optional covariates, cluster variable, bootstrap, and a generated R script that reproduces every number. |
| **Bacon, Pre-trend and HonestDiD** | Goodman-Bacon decomposition, TWFE weights, joint pre-trend tests, HonestDiD (relative magnitudes; slow, off by default, rows stream in). |
| **Continuous Treatment DiD** | Dose-response ATT(d) and marginal effect ACRT(d) with overall averages (`contdid`). |
| **Triple Differences (DDD)** | `triplediff` estimates next to ordinary DiD on eligible-only and ineligible-only units. |
| **Treatment Switching On and Off** | Dynamic effects and placebos for a binary treatment that can switch on and off (de Chaisemartin and D'Haultfoeuille, 2024); own implementation, equal to `DIDmultiplegtDYN`. |

### Estimators in "Compare DiD Estimators"

| Estimator | Implementation | Comparison group | Notes |
|---|---|---|---|
| TWFE | `fixest` | all (implicit weights) | static and event study; can be biased with staggered timing and heterogeneous effects |
| Callaway & Sant'Anna | `did` | never or not-yet treated | doubly robust / IPW / regression with covariates; analytic or multiplier-bootstrap SE |
| Sun & Abraham | own code on `fixest` | never treated | matches `fixest::sunab` exactly |
| Borusyak, Jaravel & Spiess | `didimputation` | untreated observations | post-treatment effects plus BJS pre-trend coefficients |
| Gardner two-stage | `did2s` | untreated observations | point estimates equal BJS, different SEs |
| Extended TWFE (Wooldridge) | own code on `fixest` | not yet treated or never treated | linear, Poisson, logit; covariates interacted with cohort; equals the package `etwfe` (checked in `dev/test_etwfe.R`) and, without covariates, the imputation estimators |
| Stacked DiD | own code on `fixest` | clean controls per cohort | Cengiz et al. style |

Treatment can be a 0/1 indicator (must be absorbing) or a first-treated-period variable. Standard errors are
clustered at the unit level, or at the level of an optional cluster variable. An optional **sampling weight** variable (constant
within units, for example baseline population) is used by TWFE, Callaway-Sant'Anna, Sun-Abraham, Gardner and stacked DiD;
Borusyak-Jaravel-Spiess is skipped with weights (`didimputation` has none).

**Inference**: TWFE, Sun-Abraham and stacked DiD use t(G-1) critical values and p-values with G clusters, exactly as `fixest` and
Stata's `reghdfe` do (with 12 clusters a normal approximation would report p = 0.0009 where t gives 0.0069); Callaway-Sant'Anna,
Borusyak-Jaravel-Spiess and Gardner keep their packages' normal inference. The overall table reports observations and clusters per
estimator, and the panel summary warns about few clusters and very small cohorts. On unbalanced panels the
Callaway-Sant'Anna, DDD and continuous-treatment estimators use only units observed in every period (the
unbalanced-panel code path in `did` can crash R, so it is not used); the fixed-effects estimators use all rows.

Every analysis has a collapsed **Notes and guide** group at the end of its options with two unticked options: *Method notes* and *Data and setup guide* (variables, data format, example dataset).

## Use from R (Rj editor, scripts)
Every analysis is also an R function, for example
```r
didfamily::didcompare(data = d, dep = "y", id = "id", time = "year", treat = "D")
```
and the "Show reproducible R code" option prints a standalone script (plain `fixest`, `did`, `didimputation`, `did2s`).

## References
Papers and R packages behind each analysis appear in the References section of the jamovi output (defined in `jamovi/00refs.yaml`).

## Data format
Long panel: one row per unit and period. Required: outcome, unit ID, integer-valued time period, and either a
0/1 treatment indicator or the first treated period (0/missing = never treated). Units treated in the first
observed period are dropped (no pre-period). Rows with missing values are dropped and reported.
Periods that start at zero or below, or are equally spaced by more than one (survey waves every 2 years), are
recoded internally to 1, 2, 3, ... (event time then counts periods) and shown with their original values;
unequally spaced periods are flagged. Text, dates and non-0/1 treatment indicators get a clear message.

Sample data ships with the module as `.omv` files with variable types and descriptions preset (IDs as ID, binary indicators as
nominal, periods and cohorts as continuous) and appears in jamovi under **Open > Data Library > Difference-in-Differences Family**; each entry says
which analysis to use and which variable goes where. The same data as CSV are in `dev/csv/` (regenerate the `.omv` files with
`dev/make_omv.R`):
- `sim_staggered.csv` — staggered adoption with heterogeneous dynamic effects (TWFE is visibly biased). `id, year, y, D, g`
- `sim_2x2.csv` — two periods, covariate-dependent trends, true effect 1.5. `id, year, y, treated, x, region`
- `sim_ddd.csv` — triple-difference design with a group-wide shock. `id, year, y, g, elig`
- `sim_dose.csv` — continuous treatment, ATT(d) = 2d + d², true ACRT(d) = 2 + 2d. `id, year, y, g, dose`
- `mpdta.csv` — minimum-wage example from the `did` package (use "First treated period" with `first.treat`).
- `castle.csv` — Castle Doctrine laws and crime in 50 US states, 2000-2010 (Cheng & Hoekstra 2013, via `bacondecomp`):
  staggered adoption, 29 never-treated states. `state, year, post, effyear, l_homicide, ..., pop2000` (`pop2000` = 2000 population, a time-constant weight)

## Plot options
Every analysis follows the plot theme chosen in jamovi's Plots tab for the base look. The module adds a "Plot style"
section: colour-blind-safe colours or grayscale (shades, point shapes and line types, for black-and-white printing);
for the event-study plots also error bars or shaded bands, and the spacing between estimators.

## Install (development)
Requires jamovi 2.6+ and the jamovi compiler (`jmvtools`):
```r
jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")   # adjust to your jamovi folder
```
This builds the module, bundles its R dependencies (~130 MB) and installs it into jamovi. A prebuilt
`.jmo` can be installed from jamovi's module library (sideload).

## Validation
`dev/test_*.R` hold ~300 checks: simulated designs with known effects, exact identities (Bacon weights,
TWFE weights, the fast Goodman-Bacon code vs `bacondecomp`, Sun-Abraham vs `fixest::sunab`, generated R code vs module output), the `mpdta` example,
the Castle Doctrine data (TWFE 0.0818, 91% of the weight on treated-vs-untreated comparisons), sampling weights against
replicated units (a weight of w equals w copies of the unit, for every weighted estimator), t(G-1) p-values against `fixest`, and
the variable types stored in the `.omv` files.

## Known limitations
- No dCDH estimator (`DIDmultiplegtDYN` needs `polars`, absent from jamovi's package snapshot), no synthetic DiD
  (`synthdid` is GitHub-only), no `fect` (its FEct estimator equals the imputation estimators already included).
- Wooldridge's extended TWFE is not listed separately: with never-treated controls it reproduces the imputation
  estimates, and with leads added it equals Sun-Abraham.
- No anticipation periods and no user-chosen reference period (event time is relative to e = -1).
- Sampling weights exist only in "Compare DiD Estimators" (not in the diagnostics, DDD, 2x2 or continuous-treatment analyses).
- Panel data only: repeated cross-sections of individuals with a state-level policy must first be collapsed to unit-period cells.
- No wild cluster bootstrap, no uniform (sup-t) confidence bands, no unit-specific trends.
- Treatment must be absorbing (staggered adoption). Covariates that never change within a unit are used by
  Callaway-Sant'Anna and DDD only; fixed-effects estimators absorb them.
- HonestDiD covers relative-magnitude restrictions only and needs consecutive pre- and post-treatment event times.
- Continuous-treatment curves can be biased near the smallest doses with low spline degree (see the notes in the output).
- Only labels, tables and options are translated (Japanese, Simplified Chinese); explanatory notes and messages are English.
- Package versions are those of jamovi's CRAN snapshot (for example `did` 2.3.0, `contdid` 0.1.0, `triplediff` 0.2.0).

## Notes for developers
- UI: `compilerMode: tame` in the `.u.yaml` files keeps the hand-written layout; `aggressive` would regenerate a flat one on every build.
- Layout: `R/core_*.R` is plain R with no jamovi dependency (preparation, estimators, diagnostics, DDD, 2x2,
  continuous treatment, panel checks, code generation, plots); `R/*.b.R` and `jamovi/*.yaml` are the jamovi front end.
- jamovi replaces `stats::as.formula` with a validator that rejects unlisted functions inside formulas
  (`sunab()`, `i()` built from strings, ...). Estimators therefore avoid string-built formulas with function
  calls. `dev/sandbox.R` emulates the validator so tests catch problems without opening jamovi.
- Tests: `dev/test_core.R`, `test_code.R`, `test_2x2.R`, `test_panel.R`, `test_dose.R`, `test_castle.R`
  (see each header for the library-path setup); `dev/edge_case.R` (one odd dataset per process, to catch R crashes) and
  `dev/perf.R` (timing: 100k rows take about 30 s for all estimators). Notes and findings: `dev/NOTES.md`.
- Translations: `jamovi/i18n/*.po`; `dev/fill_po.py` holds the translation table.
