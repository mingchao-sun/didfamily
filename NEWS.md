# didfamily 1.0.0

* First public release (jamovi library submission). Contains everything listed below; the earlier 0.x entries are development history.

# didfamily 0.8.0

- **New analysis: Extended TWFE (Wooldridge)** (menu DiD > Estimators). Cohort-by-period treatment effects in a regression with cohort and period effects, with covariates demeaned within cohort and interacted with the effects, for **linear, Poisson (counts) and logit (0/1) outcomes**, with not-yet-treated or never-treated controls (the latter adds pre-treatment effects); overall effect, event study and effects by cohort and calendar time, with delta-method standard errors. It is a native implementation (the R package `etwfe` cannot run inside jamovi's formula sandbox); coefficients and all four aggregations with their standard errors equal `etwfe` (McDermott, 2023) to 5e-8 in 40 comparisons (`dev/test_etwfe.R`), and with a linear model and no covariates the point estimates equal the imputation estimators. Sample data: Simulated staggered counts.

# didfamily 0.7.0

- **Pre-trend power** in the analysis *Bacon, Pre-trend and HonestDiD* (section Pre-trend power): the power of the conventional pre-trend test against a linear violation of parallel trends (Roth, 2022), the slopes the test detects with 50 and 80 percent power, the bias they would put in each post-treatment effect, and a plot. Base-R engine (`R/core_pretrend.R`, `mvtnorm` for the box probability), checked against closed forms, Monte Carlo and the rejection rate of the test in simulated panels (`dev/test_pretrend.R`).
- `didcs` and the other estimator analyses are unchanged. A note for users looking for Wooldridge's extended TWFE: it is not offered as a separate analysis because the R package `etwfe` cannot run inside jamovi (its formulas use `i()`, blocked by the formula sandbox); with never-treated controls and no covariates it equals the imputation estimator (Borusyak, Jaravel and Spiess), which is in the module.

# didfamily 0.6.0

- **New analysis: Treatment Switching On and Off** (menu DiD > Diagnostics and extensions). Dynamic effects and placebo tests of a binary treatment that can turn on and off (de Chaisemartin and D'Haultfoeuille, 2024), which the other analyses cannot handle because they need an absorbing treatment. Base-R engine (`R/core_dcdh.R`), no extra package dependencies; effects, standard errors, the number of switchers and the joint placebo test equal the R package DIDmultiplegtDYN on balanced and unbalanced panels, with groups switching on, off or both (`dev/test_dcdh.R`). Sample data: Simulated treatment that switches on and off.
- Binary treatment, no covariates or weights in this version.

# didfamily 0.5.0 (renamed)

- The module is now called **didfamily** (title "Difference-in-Differences Family"), like scmfamily, rdfamily, ivfamily and matchfamily; it was **didcompare**. The analysis "Compare DiD Estimators" keeps the function name `didcompare()`, so reproducible code now starts `didfamily::didcompare(...)`. Uninstall the old module "DiD Compare" if it is still listed; analyses saved with it need the new module. Menu group, analyses and numbers are unchanged.

# didcompare 0.5.0

- **Method notes and Data and setup guide** now follow the rdfamily / scmfamily layout: headed sections, every option on the screen explained, how to read the results, limitations; the guide gives the sample dataset, an exact field-settings table, expected results and a data-layout example. Run-specific remarks appear under "This run".
- **One analysis per estimator**: the DiD menu now has an "Estimators" sub-menu with Two-Way Fixed Effects (TWFE), Callaway & Sant'Anna,
  Sun & Abraham, Borusyak, Jaravel & Spiess (Imputation), Gardner Two-Stage (did2s) and Stacked DiD, each with its own
  options (Callaway & Sant'Anna: control group, estimation method, bootstrap, group-time heat map, cohort and calendar-time effects),
  overall ATT table, dynamic effects table with confidence intervals and p-values, event-study plot and reproducible R code.
  They run the same code as "Compare DiD Estimators" (dev/test_single.R checks the numbers are identical).
- Menu sub-headings: Estimators, Compare, Diagnostics and extensions.
- **Cleaner option panels**: every option group is collapsed by default, and every analysis ends with a collapsed "Notes and guide" group with two unticked options, "Method notes" (what the method does, how to read the results) and "Data and setup guide" (which variable goes where, data format, example dataset), as in the rdfamily and scmfamily modules. The long notes are no longer shown by default.

# didcompare 0.4.1

- **Menu**: analyses are now in workflow order (panel check, 2x2, compare, diagnostics, continuous treatment, DDD); renamed "DiD Panel Diagnostics" to "Panel Structure Check" and "DiD Diagnostics and Sensitivity" to "Bacon, Pre-trend and HonestDiD".
- **Lighter, safer build**: the Gardner two-stage estimator (from `did2s`) and the HonestDiD relative-magnitudes code are now
  vendored (MIT, see NOTICE), so the module no longer depends on `did2s`/`HonestDiD` and therefore not on CVXR, gmp, clarabel
  (Rust), highs (CMake), osqp, scs or ECOSolveR. Results are identical to the packages (dev/test_did2s.R, dev/test_honestdid.R).
- **License**: module code is now `GPL-2 | GPL-3`, so it can be combined with GPL-2-only dependencies (did, BMisc, slam, staggered)
  and GPL-3 ones (fixest). NOTICE lists the vendored code.

# didcompare 0.4.0

Result of a full review against jamovi's module requirements and economists' expectations.

- **Sampling weights** in Compare DiD Estimators (TWFE, Callaway-Sant'Anna, Sun-Abraham, Gardner, stacked DiD); verified by
  frequency-weight equivalence (weight w = w copies of the unit) and against fixest; generated R code includes weights.
- **t(G-1) inference** for TWFE, Sun-Abraham and stacked DiD, matching fixest and Stata's reghdfe; the overall table now shows
  observations and clusters per estimator; warnings for few clusters and tiny cohorts.
- **References**: papers and R packages are cited in the output (jamovi/00refs.yaml).
- **R interface**: optional variables have defaults, so each analysis can be called as an R function; analysis and option
  descriptions for the R help.
- **Example data** are now `.omv` files with variable types and descriptions (new `pop2000` weight in the Castle data); CSVs moved to
  dev/csv so they are no longer installed as R datasets (LazyData removed).
- LICENSE file (GPL-3) added.
- The CS bootstrap restores the caller's random-number state; infinite outcomes are dropped like missing values.

# didcompare 0.3.0

- Data library: six example datasets (new: Castle Doctrine laws) are installed with the module and listed under
  Open > Data Library > didcompare, each with a note on which analysis and variables to use.
- Plot style options in every analysis: colour-blind-safe colours or grayscale (shades, shapes, line types);
  error bars or shaded bands and spacing for the event-study plots.
- Fixed: variables that jamovi classes as ID (for example text state or firm names) could not be dragged into the
  unit ID box; all analyses now accept ID-type variables.
- The forest plot no longer uses a function deprecated in recent ggplot2.

# didcompare 0.2.1

- Fixed: periods starting at 0 or below were misread (a unit first treated in period 0 counted as never treated);
  equally spaced periods (every 2 years, ...) broke Sun-Abraham and stacked DiD. Periods are now recoded internally
  and shown with their original values; irregular spacing is flagged.
- Clear messages for text or date variables, a non-0/1 treatment indicator, no treated units, a constant outcome
  (also in the 2x2 analysis).
- Goodman-Bacon decomposition computed from cohort means: identical to `bacondecomp`, 0.3 s instead of minutes on
  100k rows. All estimators together take about 30 s on a 5,000-unit, 20-period panel.
- Generated R code now reseeds before each bootstrap aggregation, so it reproduces bootstrap results exactly, and
  reproduces the period recoding.
- Interface: options are grouped into collapsible sections, and dependent fields (bootstrap draws and seed, HonestDiD
  settings, Callaway-Sant'Anna extras) are disabled until they apply.

# didcompare 0.2.0

- New analyses: **DiD Panel Diagnostics** (checks, cohort table, treatment-status plot), **Classic DiD (2x2)**
  (regression, doubly robust, IPW, outcome regression), **Continuous Treatment DiD** (`contdid`).
- Compare DiD Estimators: reproducible R code output (verified to reproduce module results), forest plot,
  Callaway-Sant'Anna group-time heat map and effects by cohort / calendar time, cluster variable,
  multiplier bootstrap with seed, BJS pre-treatment coefficients.
- Callaway-Sant'Anna, DDD and continuous-treatment estimators use units observed in every period on unbalanced
  panels (did's unbalanced-panel path can crash R); the panel summary says so.
- Validation on the Castle Doctrine data; ~200 automated checks.
- Fixed: BJS failed when all covariates were time-invariant; the panel check mishandled a missing outcome variable.

# didcompare 0.1.0

First release.

- **Compare DiD Estimators**: TWFE, Callaway & Sant'Anna, Sun & Abraham, Borusyak-Jaravel-Spiess,
  Gardner two-stage and stacked DiD on one panel; overall ATT table, dynamic effects, overlaid
  event-study plot, optional covariates.
- **DiD Diagnostics and Sensitivity**: Goodman-Bacon decomposition, TWFE weights, joint pre-trend
  tests, HonestDiD (relative magnitudes).
- **Triple Differences (DDD)**: DDD next to eligible-only and ineligible-only DiD.
- Japanese and Simplified Chinese translations of labels and tables.
- Sample data: `sim_staggered.csv`, `sim_ddd.csv`, `mpdta.csv`.

Not included: dCDH (needs `polars`), synthetic DiD, continuous-treatment DiD.
