# M0 notes (2026-11-14)

Environment: jamovi 28.17.0, bundled R 4.6.0, jmvtools 28.17, Rtools45 (R 4.6 expects Rtools46, so source builds can fail).

jmvtools needs the jamovi path explicitly:
`jmvtools::prepare(home = "C:/Program Files/jamovi 28.17.0.0")` (same for `create`, `addAnalysis`, `install`).

Install compiled deps as binaries: `options(install.packages.compile.from.source = "never")`.

Feasibility (dev/check_deps.R, dev/check_deps.log):
- OK: fixest 0.14.2, did 2.5.1 (needs binary fastglm 0.1.3), DRDID, DIDmultiplegtDYN 2.4.0, didimputation, did2s, etwfe, bacondecomp, TwoWayFEWeights, fect, HonestDiD 0.2.8, contdid 0.1.2, triplediff 0.2.4
- Missing: pretrends, synthdid (GitHub-only; vendor or reimplement)
- Not yet tested: whether these deps get bundled when building the .jmo for distribution (needs a per-platform build check).

## .jmo install test (2026-11-14)
- `jmvtools::install(home=...)` works. Module installs to %APPDATA%/jamovi/modules/didfamily.
- Declaring the 10 engines in DESCRIPTION Imports bundles ~80 pkgs (~159 MB) into the module's R/ folder.
- The compiler skips deps already provided by jamovi's base and jmv module libs (ggplot2, cowplot, Formula, plyr, abind, Rcpp, data.table, cli...). Loading with R_LIBS=<module>;<base>;<jmv> works: all 10 engines + didcompare load, and did::att_gt on mpdta gives simple ATT -0.0400.
- Loading with base + module only fails (Formula/plyr/cowplot/abind missing), so the module relies on jmv being installed (it ships with jamovi).
- Not tested: opening the analysis in the jamovi GUI; Mac/Linux builds.

## GUI check (2026-11-14)
- jamovi 28.17 shows the module as "didcompare" in the Analyses ribbon with "DiD Data Diagnostics"; the analysis opens and runs (stub still has the template UI: dep/group/alt/varEq).
- In-app .libPaths(): module R, Resources/modules/jmv/R, base/R, Frameworks/R/library (matches the compiler's assumption).
- Bundled versions come from the CRAN snapshot jamovi uses (2026-05-11), not the latest: did 2.3.0 (not 2.5.1), fixest 0.14.1, DRDID 1.2.3, DIDmultiplegtDYN 2.3.3, didimputation 0.5.1, did2s 1.2.1, etwfe 0.6.2, bacondecomp 0.1.1, HonestDiD 0.2.8, fect 2.4.1. Code must target these API versions.

## M1 (2026-11-14): core layer + comparison analysis
Files: R/core_prepare.R, R/core_estimators.R (TWFE, CS, SA), R/core_plot.R, R/didcompare.b.R; jamovi/didcompare.{a,r,u}.yaml; data/sim_staggered.csv; tests: dev/test_core.R.
Run tests in a jamovi-like env (R_LIBS_USER=notthere; R_LIBS=<module R>;<jamovi>/Resources/modules/jmv/R;<jamovi>/Resources/modules/base/R).

Findings that constrain later milestones:
- jamovi replaces stats::as.formula with jmvcore::asFormula, which rejects formula function calls outside a whitelist ("Security violation: Unsafe function"). fixest::sunab() is blocked inside jamovi (works in plain R). i() interactions pass. Do not bypass; implement estimators with plain i()/data columns instead. Expect the same check for other packages that build formulas with as.formula (watch DIDmultiplegtDYN, did2s, didimputation, etwfe, fect in M2).
- Own SA implementation (i(rel, factor(cohort), ref=-1, ref2="0") + cohort-size weights) equals fixest::sunab to 1e-6 incl. SE.
- jamovi's ggtheme object carries its own colour scale; add custom scales after `+ ggtheme`.
- Dropbox can leave build/R*/00LOCK-temp behind: delete it if jmvtools::install fails with "failed to lock directory".
- Delete jamovi/0000.yaml analysis entries by hand when removing an analysis (stale entries break the namespace export).
- Sim check (true effect on treated unit-periods 1.685): TWFE 1.446 (pre-trend -0.31 at e=-4), CS 1.640, SA 1.640.

## M2 (2026-11-14): BJS, Gardner, stacked DiD
Added dc_bjs, dc_gardner, dc_stacked (R/core_estimators.R), UI/tables/notes, dev/sandbox.R (emulates jamovi's as.formula validator in plain R; test_core.R runs every estimator under it and requires results equal to plain R).
Sim check (true 1.685): TWFE 1.446, CS 1.640, SA 1.640, BJS 1.650, Gardner 1.650, Stacked 1.588 (in-memory sim); GUI CSV data: BJS = Gardner = 1.607.

Decisions / blockers:
- didimputation::did_imputation(horizon=) and did2s second_stage=~i(...) build i() formulas -> blocked by jamovi's validator. Workarounds: BJS post horizons via wname weights (weight 1 on untreated rows, indicator of e==h on treated), Gardner via hand-built event-time dummy columns. BJS pre-period coefficients are not reported (BJS use a separate pre-trend regression).
- etwfe is blocked the same way, and redundant: verified ETWFE post-cells with never-treated controls reproduce the BJS overall ATT exactly (1.606525); with leads it equals Sun-Abraham. Not listed as a separate estimator.
- DIDmultiplegtDYN (dCDH) 2.3.3 needs the 'polars' package (r-universe only, not in jamovi's CRAN snapshot) -> not usable. Options for later: own implementation of DID_l for absorbing binary treatment, or wait for a CRAN polars.
- Do not write "\1" in bash heredocs here: the backslash gets collapsed; use the Edit tool for R code with regex backreferences.

## M3 (2026-11-14): diagnostics and sensitivity
New analysis `didiag` ("DiD Diagnostics and Sensitivity"): R/core_diag.R (dc_bacon, dc_twfe_weights, dc_pretests, dc_honest_prepare/one), R/didiag.b.R, jamovi/didiag.{a,r,u}.yaml, plots in R/core_plot.R. dc_sa refactored into dc_sa_core (event-time estimates + full covariance). dc_info_html shared by both analyses. Tests: dev/test_core.R (all under the emulated sandbox).
- bacondecomp::bacon works under the sandbox (balanced panels only; implied TWFE = TWFE coefficient to 1e-6).
- TWFE weights implemented directly (residual of D on unit+period FE, normalised over treated cells); exact identity test: beta_TWFE = sum w_it tau_it.
- Pre-trend Wald tests: TWFE and SA from coefficient covariance; CS from did's att_gt W/Wpval (df = pre cells with group > t+1, consistent with p). CS test ignores the event window.
- HonestDiD 0.2.8 works under the sandbox but is slow (~8 s per Mbar value). Run one Mbar at a time with private$.checkpoint(): rows stream into the table and the run can be cancelled. Off by default. Relative-magnitudes only; smoothness (DeltaSD, ~20 s for 4 values) and pretrends power (package GitHub-only) not included. Needs consecutive pre (e < -1) and post (e >= 0) event times. CS covariance comes from did's influence functions (crossprod(IF)/n^2).
- YAML gotcha: a colon inside an unquoted title breaks the compiler ("bad indentation of a mapping entry").
- The jamovi compiler rewrites .u.yaml into a flat layout when stage: 0; nested layouts (e.g. enable: under a CheckBox) are not preserved.

## M4 (2026-11-14): covariates, DDD, translations, docs, packaging
- Covariates: dc_prepare(covs=) (numeric as is, factors -> 0/1 dummies, labels kept in info$covs, `varying` flag). Formulas are built with eval(bquote(...)) from templates + internal names dc_x1.. (string-built formulas with i() are blocked by the sandbox). Time-invariant covariates are dropped for fixed-effects estimators (collinear with unit FE; BJS crashed on them) and only CS/DDD use them (dc_cov_expr(d, fe = FALSE)). CS gets est_method dr/ipw/reg. Test: covariate-dependent trends, CS 1.230 -> 1.041 (true 1.0). Diagnostics analysis ignores covariates.
- DDD: R/core_ddd.R, analysis `dddcompare`. triplediff in jamovi's snapshot is 0.2.0 (CRAN has 0.2.4); API used: ddd(..., base_period = "universal"), agg_ddd(type = "simple"/"eventstudy") -> $aggte_ddd. Compared with CS DiD on eligible-only and ineligible-only units. sim_ddd: DDD 1.511 (true 1.489), DiD eligible 2.319, DiD ineligible 0.804 (= group shock 0.8). Partition variable must be 0/1 and constant within unit.
- Dependencies trimmed to what is used: fixest, did, didimputation, did2s, bacondecomp, HonestDiD, triplediff (+ggplot2). Build ~128 MB installed, 69 MB .jmo (didcompare_0.1.0.jmo in the project root). Clean rebuild: delete build/R4.6.0-x64-win64 first (Dropbox may keep a folder busy; a partial rm still worked).
- i18n: jamovi/i18n/{ja,zh_Hans}.po (87 strings, filled by dev/fill_po.py). Only yaml labels are extracted; R strings (notes, info, errors) are not wrapped for translation. Not verified in the GUI (would require switching jamovi's language).
- Docs: README.md, README.zh.md; sample data in data/ (sim_staggered, sim_ddd, mpdta).
- Not done / decisions pending: contdid and synthdid (not probed; synthdid GitHub-only), dCDH (polars), (author/maintainer/license and git repository: done in the release step), .omv sample files, i18n for explanatory notes, jamovi library submission, Mac/Linux build check.
- Tests: dev/test_core.R, 96 checks, all pass (log in dev/test_core.log).

## Round 5 (2026-11-14): closing the gaps found in the plan audit
Added: Classic DiD (2x2) [core_2x2.R, did22], DiD Panel Diagnostics [core_panel.R, didsetup], Continuous Treatment DiD [core_dose.R, dosedid, contdid], reproducible R code [core_code.R], forest plot + group-time heat map + cohort/calendar tables, cluster variable, CS bootstrap (reps/seed), BJS pre-period coefficients (lead dummies on untreated obs), Castle Doctrine validation.
- did 2.3.0 with allow_unbalanced_panel = TRUE can segfault R (heap corruption, input-dependent, not deterministic). Never use it: dc_balanced_subset() keeps units observed in every period; the info panel says how many were set aside.
- DRDID repeated cross-section mode needs an id column even when unused (pass a row id). contdid in the jamovi snapshot is 0.1.0 (CRAN has 0.1.2); API used: cont_did(..., target_parameter = "level", aggregation = "dose", dose_est_method = "parametric").
- Parametric dose-response is biased near the smallest doses (controls at dose 0); documented in the output notes.
- Generated R script (core_code.R) reproduces TWFE, CS, SA, BJS, Gardner and stacked ATT/SE exactly (dev/test_code.R, runs in plain R; it uses sunab()/i() which the jamovi sandbox blocks).
- Castle Doctrine: TWFE 0.0818 (se 0.0589); Bacon weights 0.908 treated-vs-untreated, 0.060 earlier-vs-later, 0.032 later-vs-earlier; CS 0.110, SA 0.110, BJS/Gardner 0.080, stacked 0.080. did cannot compute its pre-test Wald statistic there (small cohorts) -> the module reports a note.
- Not done: anticipation, user-chosen reference period, plot style options, HonestDiD smoothness, dCDH (polars), synthdid (GitHub only), fect/pretrends, Stata cross-check, .omv files, testthat layout, translating notes/messages, Mac/Linux build, library submission.
- Dropbox sometimes locks build/js/*.src.js (EPERM) during jmvtools::install: just retry.
- Heredocs with backslash sequences get mangled in the Bash tool: write Python patch scripts to files instead.

## Round 6 (2026-11-14): detail review
- dev/edge_case.R runs ~20 odd datasets, one process each (no R crash found; found real bugs): time <= 0 made first-treated-at-0 units "never treated" (0 is the internal never code) and period gaps broke SA/stacked -> dc_prepare recodes t to 1..T when t0 <= 0 or equal step > 1 (info$time_map; dc_orig_time() for display; generated code mirrors it). Friendly stops for text/date/constant/no-treated inputs.
- jamovi `compilerMode: aggressive` in .u.yaml regenerates a flat layout on every build and drops `enable:`; `tame` keeps hand-written layouts.
- Performance (dev/perf.R, 5,000 units x 20 periods): TWFE 1.4 s, CS 4.3, SA 4.3, BJS 9.5, Gardner 2.9, stacked 1.8, weights 1.0, pretests 7.7. bacondecomp::bacon took minutes (one lm with unit dummies per 2x2): replaced by dc_bacon_2x2 (cohort means; weights use unit-periods inside the window as group sizes, as bacondecomp does) -> 0.3 s, identical output (4 designs).
- Generated CS code must reseed before every aggte() (the module does).

## Round 7 (2026-11-14): data library and plot style
- Data library = `datasets:` section in jamovi/0000.yaml (name, path, description) + files in data/; jmc copies them into the module. Descriptions are not translated (only names are extracted).
- Plot style options (plotPalette, plotCI, plotDodge) are in the clearWith of every result item: a style change is a clean re-run. jamovi's own theme (Plots tab) arrives as ggtheme; our plots add it first and their colour scales last.
- BUG found via the Castle data: ID-type variables (text with many levels, key icon) were greyed out and undroppable because the VariableSupplier in each .u.yaml lacked `permitted: [numeric, factor, id]` (jmv sets it on the supplier; the option-level `permitted` is not enough). Fixed in all six analyses.
- geom_errorbarh is deprecated in ggplot2 4.x: use geom_errorbar(orientation = "y").
- dev/test_plots.R renders all plots in both palettes and CI styles (no warnings allowed).

## Round 8 (2026-11-15): full review (jamovi requirements, economists' standards)
Fixed: no refs/descriptions/LICENSE; optional variables lacked `default: null` (R-function call failed with "argument covs is missing");
CSVs in data/ were installed as R datasets (castle, mpdta ... masked did::mpdta; LazyData); normal instead of t(G-1) inference for fixest-based
estimators (12 clusters: p 0.0009 vs fixest 0.0069); no N in the results table; no sampling weights; set.seed leaked into the caller's RNG;
infinite outcomes crashed the new constant-outcome check; generated code broke when a weights column had a trailing comment.
Open (need a decision or are out of reach here):
- License: did, BMisc, slam, staggered are GPL-2 only (R's "GPL-2" = version 2 only); fixest etc. are GPL-3: no single GPL version satisfies both when redistributed together. Options: "GPL-2 | GPL-3" for the module's own code, or ask the jamovi library maintainers.
- Cross-platform build risk: HonestDiD -> CVXR -> gmp, Rglpk (GLPK), clarabel (Rust toolchain), highs (CMake), osqp, scs: 38 packages with compiled code, 72 MB. Only Windows verified. If the jamovi library's Mac/Linux builds lack these, the whole module fails; fallback = make HonestDiD optional or split it into a second module.
- Not implemented: repeated cross-section inputs, wild cluster bootstrap, uniform confidence bands, unit-specific trends, anticipation, user-chosen reference period, dCDH (polars), synthetic DiD, weights outside the comparison analysis, generated code outside the comparison analysis.
- Docs are English only for notes/messages; dataset descriptions deliberately untranslated.

## Round 9: license and build-risk decisions (user: "decide now")
- Vendored did2s two-stage estimator (R/core_did2s.R) and the HonestDiD RM subset (R/vendor_honestdid.R, generated by dev/gen_honestdid.R);
  equality tests vs the packages pass. Imports no longer contain did2s/HonestDiD; CVXR chain gone. Rglpk (GLPK) and lpSolveAPI remain
  (HonestDiD's conditional LP); both are small and self-contained. Generated user scripts still call did2s (user installs it).
- License set to "GPL-2 | GPL-3" (NOTICE). Revert in DESCRIPTION if a single GPL-3 is preferred, but then did/BMisc/slam/staggered conflict.
## Ribbon icon
jamovi's client bundles ribbon icons per module name (CSS rule .jmv-ribbon-button[data-name=menuGroup], lower case; ours is "did"), and a module .jmo cannot carry its own. icon/did.svg (45x30, same size as the built-in ones) is the proposed icon for the DiD menu group; it has to be supplied to the jamovi team when the module is submitted to the library.

## Round 10: single-estimator analyses (0.5.0)
- Six new analyses (didtwfe, didcs, didsa, didbjs, didgardner, didstacked) in the menu sub-group "Estimators", asked for after comparing with scm / rdd, which list each method plus a comparison entry.
- They are generated: `python dev/gen_single.py` writes jamovi/did<key>.{a,u,r}.yaml and R/did<key>.b.R; the logic is R/core_single.R (dc_single_run / dc_single_plot), which calls the dc_* estimators of core_estimators.R and dc_code_compare for the R script. Edit the generator or core_single.R, never the outputs, and re-run the generator.
- BJS has no weights option (didimputation has none). Only Callaway & Sant'Anna has estimator-specific options so far; ideas for later: anticipation, reference period, wild bootstrap, uniform bands.
- Translations: `jmvtools::i18nUpdate(home=...)` then `python dev/fill_po.py`.
- Option panels: all CollapseBox collapsed by default; last box "Notes and guide" with showHelp (notes item, visible only when ticked) and showGuide (guide item, text from R/core_guide.R, one entry per analysis). dev/add_notes_guide.py was the one-off migration of the six hand-written analyses (do not run again); the generated ones get the same from dev/gen_single.py. The panel check's empty-state message moved from `notes` to its own `info` item.
- Fix after first use: the notes were only filled at the end of a successful run and the guide sat below the plots, so ticking the options before assigning variables showed nothing (the text was far down, behind the empty plot area). Now the guide item is first in every results file, and the notes are set at the start of every run from static text (dc_method_html in R/core_guide.R, dc_single_notes for the estimators), then replaced by the full notes after a real run.
