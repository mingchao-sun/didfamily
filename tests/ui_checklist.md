# didfamily: checklist for testing in jamovi (new in 0.6.0)

Menu: DiD. Install: `jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`.
Open sample data from File > Open > Data Library > Difference-in-Differences Family.

## 1. Treatment Switching On and Off (DiD > Diagnostics and extensions)
Dataset: Simulated treatment that switches on and off. Outcome `y`; Group ID `id`; Time period `year`; Treatment (0/1) `D`; Number of effects 3; Number of placebos 2.
Set `year` to Continuous and `id` to ID or Nominal; `D` Continuous or Nominal (both work).
Expect: a note "300 groups and 10 periods. 218 groups change their treatment (218 switch on from 0, 0 switch off from 1); 82 never change it."
Effects: Effect_1 1.995 (s.e. 0.107), Effect_2 1.891 (0.119), Effect_3 1.168 (0.134), 218 switchers each; Placebo_1 -0.096 (0.112), Placebo_2 0.026 (0.113); the "Periods from the switch" column reads 0, 1, 2, -1, -2.
Joint tests: placebos p = 0.57; effects p < .001. The plot shows placebos (orange triangles, negative x) and effects (blue circles) with a dashed zero line.
Then: effects 5 and placebos 3; confidence level 90; Plot colors "Grayscale"; untick Event study plot; tick Method notes and Data and setup guide.
Messages to try: a treatment with values 0/1/2 ("must be binary"); a time variable with decimals; no variable switching (set a copy of D to all 0).

## 2. The older analyses (spot check after the 0.6.0 reinstall)
- Compare DiD Estimators on "Minimum wage and teen employment (mpdta)": Outcome `lemp`, Unit `countyreal`, Time `year`, Treatment `first.treat` as First treated period; the table has six estimators; TWFE -0.0365.
- Panel Structure Check on the same data lists the cohorts 2004, 2006, 2007.
- One estimator from the Estimators sub-menu (for example Callaway & Sant'Anna) gives overall ATT -0.040.

## Extended TWFE (Wooldridge) (DiD > Estimators)
Dataset: Simulated staggered adoption. Outcome `y`; Unit ID `id`; Time `year`; Treatment variable `g` with "First treated period"; defaults (linear, not yet treated).
Expect: overall ATT 1.607 (s.e. 0.024) - the same as Borusyak-Jaravel-Spiess; the event study 0.559 at e = 0, 1.139, 1.734; the plot as for the other estimators. Control group "Never treated": ATT 1.637 (s.e. 0.038) and pre-treatment effects 0.066 and 0.044 at e = -3 and -2 (the pre-trend check). Tick "Effects by cohort" and "Effects by calendar time": two tables.
Dataset: Simulated staggered counts. Outcome `y`; Unit ID `id`; Time `year`; Treatment `g` (first treated period); Covariate `x1`; Model for the outcome "Poisson"; Event window from -3 to 4.
Expect: overall effect 2.601 (s.e. 0.144) additional events; by event time 1.759, 2.224, 2.777 at e = 0, 1, 2; by cohort 2.862, 2.539, 2.140. Never treated: 2.506 (s.e. 0.184) and leads -0.124, 0.062. Linear model: 2.722 (s.e. 0.130). Show reproducible R code: a script with `etwfe(...)` and four `emfx(...)` calls (it needs the package etwfe in ordinary R). Method notes and guide on/off.
Messages: model Logit with a non-0/1 outcome -> a note in the table; "Never treated" on data without never-treated units -> a note; Poisson with a negative outcome -> a note.

## Pre-trend power (Bacon, Pre-trend and HonestDiD > Pre-trend power)
Dataset: Simulated staggered adoption. Outcome `y`; Unit ID `id`; Time `year`; Treatment `g` (first treated period); untick Bacon, weights and pre-trend tests if you want only this; tick "Power of the pre-trend test"; event study Sun & Abraham; hypothesized trend 0.1.
Expect: largest |t| 1.598 (the test passes); joint Wald p 0.362; size 0.123; power against a trend of 0.1: 1.000; slope detected with 50% power 0.026 and with 80% power 0.040. The second table lists e = 0 to 4 (estimates 0.611, 1.157, 1.774, 2.311, 1.954) with the bias of the 50% slope (0.026, 0.053, 0.079, 0.106, 0.132), of the 80% slope and of the hypothesized slope (0.1, 0.2, ... 0.5). The plot shows the event study with three straight lines through e = -1.
Then: event study Callaway & Sant'Anna (almost the same numbers: slopes 0.026 and 0.040); hypothesized trend 0 (the "Power against a trend" row and the last column disappear); event window from -2 (one pre-period: size 0.05); from -1 (a message that a pre-period is needed). Method notes on/off.

## Things to look at in the interface
- The new analysis sits under DiD > Diagnostics and extensions with the title "Treatment Switching On and Off"; collapsed groups open and close.
- Tables fill in as variables are dropped into the boxes and clear when one is removed; the hint text appears before all four variables are chosen.
- Japanese and Chinese interface (the module ships ja and zh translations): the new analysis shows English text until translations are added to `jamovi/i18n/*.po` (expected).
- Number formats, p-values, the plot with the dark theme and another jamovi palette.
