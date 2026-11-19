# Generates the six single-estimator analyses (a/u/r yaml + thin .b.R) for the DiD menu.
# Run from the project root:  python dev/gen_single.py
# The shared logic lives in R/core_single.R; edit that file or this generator, not the outputs.
import os, textwrap

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

SPECS = [
    dict(name='didtwfe', key='twfe', title='Two-Way Fixed Effects (TWFE)', short='TWFE',
         desc="Two-way fixed effects regression of the outcome on the treatment indicator with unit and period fixed effects, static and event study, with a warning that TWFE can be biased under staggered adoption.",
         weights=True, cs=False),
    dict(name='didcs', key='cs', title="Callaway & Sant'Anna", short="Callaway-Sant'Anna",
         desc="Callaway and Sant'Anna group-time effects with never-treated or not-yet-treated controls, doubly robust, IPW or outcome-regression estimation, multiplier bootstrap, and aggregation to overall, event-time, cohort and calendar effects.",
         weights=True, cs=True),
    dict(name='didsa', key='sa', title='Sun & Abraham', short='Sun-Abraham',
         desc="Sun and Abraham interaction-weighted estimator with cohort-by-event-time effects and never-treated controls, with overall and event-study effects.",
         weights=True, cs=False),
    dict(name='didbjs', key='bjs', title='Borusyak, Jaravel & Spiess (Imputation)', short='Borusyak-Jaravel-Spiess',
         desc="Borusyak, Jaravel and Spiess imputation estimator: unit and period effects are fitted on untreated observations and used to impute treated counterfactuals, with a pre-trend test on untreated observations.",
         weights=False, cs=False),
    dict(name='didgardner', key='gardner', title='Gardner Two-Stage (did2s)', short='Gardner (did2s)',
         desc="Gardner two-stage difference-in-differences: untreated outcomes are residualised on unit and period effects, then regressed on treatment, with corrected clustered standard errors.",
         weights=True, cs=False),
    dict(name='didetwfe', key='etwfe', title='Extended TWFE (Wooldridge)', short='Extended TWFE (Wooldridge)',
         desc="Wooldridge's extended two-way fixed effects: cohort-by-period treatment effects with covariates interacted with cohort, for linear, Poisson and logit models, aggregated to an overall effect, an event study and effects by cohort and calendar time.",
         weights=False, cs=False, et=True),
    dict(name='didstacked', key='stacked', title='Stacked DiD', short='Stacked DiD',
         desc="Stacked difference-in-differences: one sub-experiment per treatment cohort with clean controls, stacked with unit-by-experiment and period-by-experiment fixed effects.",
         weights=True, cs=False),
]
for _sp in SPECS:
    _sp.setdefault('et', False)

def opt(name, title, typ, default=None, desc='', extra=''):
    s = f"    - name: {name}\n      title: {title}\n      type: {typ}\n"
    if extra: s += extra
    if default is not None: s += f"      default: {default}\n"
    s += f"      description:\n          R: {desc}\n\n"
    return s

def a_yaml(sp):
    n, t = sp['name'], sp['title']
    s = f"---\nname:  {n}\ntitle: {t}\nmenuGroup: DiD\nmenuSubgroup: Estimators\nmenuTitle: {t}\nversion: '1.0.0'\njas: '1.2'\n\n"
    s += "description:\n    main: '" + sp['desc'].replace("'", "''") + "'\n    R:\n        dontrun: true\n        usage: |\n"
    s += "            # d is a data frame with one row per unit and period\n"
    s += f"            {n}(data = d, dep = 'y', id = 'id', time = 'year', treat = 'D')\n\n"
    s += "options:\n    - name: data\n      type: Data\n\n"
    s += opt('dep', 'Outcome', 'Variable', None, "'a string naming the outcome variable in `data`'", "      permitted: [numeric]\n")
    s += opt('id', 'Unit ID', 'Variable', None, "'a string naming the unit identifier variable in `data`'", "      permitted: [numeric, factor, id]\n")
    s += opt('time', 'Time period', 'Variable', None, "'a string naming the period variable in `data`; values must be integers (a year or a period number)'", "      permitted: [numeric]\n")
    s += opt('treat', 'Treatment variable', 'Variable', None, "'a string naming the treatment variable in `data`: a 0/1 indicator that never turns off, or the first treated period (see `treatType`)'", "      permitted: [numeric, factor, id]\n")
    s += opt('covs', 'Covariates', 'Variables', 'null', "'a vector of strings naming covariates in `data` (optional); text and factor variables are dummy-coded'", "      permitted: [numeric, factor]\n")
    s += opt('cluster', 'Cluster variable (optional)', 'Variable', 'null', "'a string naming a cluster variable that is constant within units (optional); by default standard errors are clustered by unit'", "      permitted: [numeric, factor, id]\n")
    if sp['weights']:
        s += opt('weights', 'Sampling weights (optional)', 'Variable', 'null', "'a string naming a positive weight variable that is constant within units (optional), for example a baseline population'", "      permitted: [numeric]\n")
    s += ("    - name: treatType\n      title: Treatment variable is\n      type: List\n      options:\n"
          "        - title: Treatment indicator D (0/1)\n          name: indicator\n"
          "        - title: First treated period (0 or missing = never)\n          name: cohort\n"
          "      default: indicator\n      description:\n"
          "          R: '`''indicator''` (default) if `treat` is a 0/1 treatment indicator, or `''cohort''` if it holds the first treated period (0 or missing = never treated)'\n\n")
    if sp['cs']:
        s += ("    - name: csControl\n      title: Control group\n      type: List\n      options:\n"
              "        - title: Never treated\n          name: nevertreated\n"
              "        - title: Not yet treated\n          name: notyettreated\n"
              "      default: nevertreated\n      description:\n"
              "          R: '`''nevertreated''` (default) or `''notyettreated''`, the control group'\n\n")
        s += ("    - name: csMethod\n      title: Estimation method (with covariates)\n      type: List\n      options:\n"
              "        - title: Doubly robust\n          name: dr\n"
              "        - title: Inverse probability weighting\n          name: ipw\n"
              "        - title: Outcome regression\n          name: reg\n"
              "      default: dr\n      description:\n"
              "          R: '`''dr''` (default), `''ipw''` or `''reg''`, the estimation method used with covariates'\n\n")
        s += opt('bootstrap', 'Multiplier bootstrap standard errors', 'Bool', 'false', "'`TRUE` or `FALSE` (default), use the multiplier bootstrap for standard errors'")
        s += opt('bootReps', 'Bootstrap draws', 'Integer', '999', "'the number of bootstrap draws (default 999)'", "      min: 99\n      max: 9999\n")
        s += opt('seed', 'Random seed', 'Integer', '1234', "'the random seed used for the bootstrap (default 1234)'", "      min: 1\n      max: 2147483647\n")
        s += opt('showGT', 'Group-time effects (heat map)', 'Bool', 'false', "'`TRUE` or `FALSE` (default), show the group-time effects as a heat map'")
        s += opt('aggGroup', 'Effects by cohort', 'Bool', 'false', "'`TRUE` or `FALSE` (default), tabulate effects by cohort'")
        s += opt('aggCalendar', 'Effects by calendar time', 'Bool', 'false', "'`TRUE` or `FALSE` (default), tabulate effects by calendar time'")
    if sp['et']:
        s += ("    - name: etControl\n      title: Control group\n      type: List\n      options:\n"
              "        - title: Not yet treated\n          name: notyettreated\n"
              "        - title: Never treated\n          name: nevertreated\n"
              "      default: notyettreated\n      description:\n"
              "          R: '`''notyettreated''` (default) uses units not yet treated as well as never-treated units, `''nevertreated''` only never-treated units (and adds pre-treatment effects)'\n\n")
        s += ("    - name: etFamily\n      title: Model for the outcome\n      type: List\n      options:\n"
              "        - title: Linear (any outcome)\n          name: gaussian\n"
              "        - title: Poisson (counts, non-negative outcomes)\n          name: poisson\n"
              "        - title: Logit (0/1 outcomes)\n          name: logit\n"
              "      default: gaussian\n      description:\n"
              "          R: '`''gaussian''` (default), `''poisson''` or `''logit''`, the model; effects are reported on the scale of the outcome'\n\n")
        s += opt('aggGroup', 'Effects by cohort', 'Bool', 'false', "'`TRUE` or `FALSE` (default), tabulate effects by cohort'")
        s += opt('aggCalendar', 'Effects by calendar time', 'Bool', 'false', "'`TRUE` or `FALSE` (default), tabulate effects by calendar time'")
    s += opt('showHelp', 'Method notes', 'Bool', 'false', "'`TRUE` or `FALSE` (default), show notes on what the method does and how to read the results'")
    s += opt('showGuide', 'Data and setup guide', 'Bool', 'false', "'`TRUE` or `FALSE` (default), show a guide to the variables and the data format'")
    s += opt('showCode', 'Show reproducible R code', 'Bool', 'false', "'`TRUE` or `FALSE` (default), show an R script that reproduces the results'")
    s += opt('minE', 'Event window, from', 'Integer', '-4', "'the most negative event time shown (default -4)'", "      min: -50\n      max: -1\n")
    s += opt('maxE', 'Event window, to', 'Integer', '4', "'the largest event time shown (default 4)'", "      min: 0\n      max: 50\n")
    s += opt('ciLevel', 'Confidence level (%)', 'Number', '95', "'the confidence level in percent (default 95)'", "      min: 50\n      max: 99.9\n")
    s += ("    - name: plotPalette\n      title: Plot colors\n      type: List\n      options:\n"
          "        - title: Colors (colour-blind safe)\n          name: color\n"
          "        - title: Grayscale (for printing)\n          name: gray\n"
          "      default: color\n      description:\n"
          "          R: '`''color''` (default, colour-blind safe) or `''gray''` (shades, shapes and line types for printing)'\n\n")
    s += ("    - name: plotCI\n      title: Show confidence intervals as\n      type: List\n      options:\n"
          "        - title: Error bars\n          name: bars\n"
          "        - title: Shaded bands\n          name: ribbon\n"
          "      default: bars\n      description:\n"
          "          R: '`''bars''` (default) or `''ribbon''`, draw confidence intervals as error bars or shaded bands'\n\n")
    s += "...\n"
    return s

def u_yaml(sp):
    n, t = sp['name'], sp['title']
    s = f"title: {t}\nname: {n}\njus: '3.0'\nstage: 0\ncompilerMode: tame\nchildren:\n"
    def box(label, name, maxone=True):
        r = f"      - type: TargetLayoutBox\n        label: {label}\n        children:\n          - type: VariablesListBox\n            name: {name}\n"
        if maxone: r += "            maxItemCount: 1\n"
        r += "            isTarget: true\n"
        return r
    s += "  - type: VariableSupplier\n    persistentItems: false\n    permitted:\n      - numeric\n      - factor\n      - id\n    stretchFactor: 1\n    children:\n"
    s += box('Outcome', 'dep') + box('Unit ID', 'id') + box('Time period', 'time') + box('Treatment variable', 'treat')
    s += box('Covariates (optional)', 'covs', False) + box('Cluster variable (optional)', 'cluster')
    if sp['weights']:
        s += box('Sampling weights (optional)', 'weights')
    s += "  - type: LayoutBox\n    margin: large\n    children:\n      - type: ComboBox\n        name: treatType\n"
    if sp['cs']:
        s += ("  - type: CollapseBox\n    label: Settings\n    collapsed: true\n    stretchFactor: 1\n    children:\n"
              "      - type: LayoutBox\n        margin: large\n        children:\n"
              "          - type: ComboBox\n            name: csControl\n"
              "          - type: ComboBox\n            name: csMethod\n"
              "          - type: CheckBox\n            name: bootstrap\n            children:\n"
              "              - type: TextBox\n                name: bootReps\n                format: number\n                enable: (bootstrap)\n"
              "              - type: TextBox\n                name: seed\n                format: number\n                enable: (bootstrap)\n")
    if sp['et']:
        s += ("  - type: CollapseBox\n    label: Settings\n    collapsed: true\n    stretchFactor: 1\n    children:\n"
              "      - type: LayoutBox\n        margin: large\n        children:\n"
              "          - type: ComboBox\n            name: etControl\n"
              "          - type: ComboBox\n            name: etFamily\n")
    s += ("  - type: CollapseBox\n    label: Event window and confidence level\n    collapsed: true\n    stretchFactor: 1\n    children:\n"
          "      - type: LayoutBox\n        margin: large\n        children:\n"
          "          - type: TextBox\n            name: minE\n            format: number\n"
          "          - type: TextBox\n            name: maxE\n            format: number\n"
          "          - type: TextBox\n            name: ciLevel\n            format: number\n")
    s += "  - type: CollapseBox\n    label: Additional output\n    collapsed: true\n    stretchFactor: 1\n    children:\n      - type: LayoutBox\n        margin: large\n        children:\n"
    if sp['cs']:
        s += "          - type: CheckBox\n            name: showGT\n          - type: CheckBox\n            name: aggGroup\n          - type: CheckBox\n            name: aggCalendar\n"
    if sp['et']:
        s += "          - type: CheckBox\n            name: aggGroup\n          - type: CheckBox\n            name: aggCalendar\n"
    s += "          - type: CheckBox\n            name: showCode\n"
    s += ("  - type: CollapseBox\n    label: Plot style\n    collapsed: true\n    stretchFactor: 1\n    children:\n"
          "      - type: LayoutBox\n        margin: large\n        children:\n"
          "          - type: ComboBox\n            name: plotPalette\n"
          "          - type: ComboBox\n            name: plotCI\n")
    s += ("  - type: CollapseBox\n    label: Notes and guide\n    collapsed: true\n    stretchFactor: 1\n    children:\n"
          "      - type: LayoutBox\n        margin: large\n        children:\n"
          "          - type: CheckBox\n            name: showHelp\n"
          "          - type: CheckBox\n            name: showGuide\n")
    return s

REFS = dict(twfe=[], cs=['callaway2021', 'did'], sa=['sun2021', 'fixest'], bjs=['borusyak2024', 'didimputation'],
            gardner=['gardner2022', 'did2s'], stacked=['cengiz2019', 'fixest'], etwfe=['wooldridge2021', 'wooldridge2023', 'etwfe'])

def r_yaml(sp):
    n, t, k = sp['name'], sp['title'], sp['key']
    cw = "plotPalette, plotCI, dep, id, time, treat, treatType, covs, cluster" + (", weights" if sp['weights'] else "")
    if sp['cs']:
        cw += ", csControl, csMethod, bootstrap, bootReps, seed"
    if sp['et']:
        cw += ", etControl, etFamily"
    s = f"---\nname:  {n}\ntitle: {t}\njrs:   '1.1'\n\nitems:\n"
    s += "    - name: guide\n      title: Data and setup guide\n      type: Html\n      visible: (showGuide)\n\n"
    s += "    - name: info\n      title: Panel summary\n      type: Html\n\n"
    s += f"    - name: summary\n      title: Overall ATT\n      type: Table\n      rows: 0\n      clearWith: [{cw}, ciLevel]\n      columns:\n"
    cols = [('method', 'Estimator', 'text', ''), ('att', 'ATT', 'number', ''), ('se', 'SE', 'number', ''),
            ('lower', 'Lower', 'number', '          superTitle: Confidence interval\n'),
            ('upper', 'Upper', 'number', '          superTitle: Confidence interval\n'),
            ('p', 'p', 'number', '          format: zto,pvalue\n'),
            ('control', 'Comparison group', 'text', ''), ('nobs', 'Observations', 'integer', ''),
            ('ncl', 'Clusters', 'integer', ''), ('note', 'Note', 'text', '')]
    for c, ti, ty, ex in cols:
        s += f"        - name: {c}\n          title: {ti}\n          type: {ty}\n{ex}"
    if REFS[k]:
        s += "      refs:\n" + "".join(f"        - {r}\n" for r in REFS[k])
    s += "\n"
    s += f"    - name: dynamic\n      title: Dynamic effects by event time (reference period = -1)\n      type: Table\n      rows: 0\n      clearWith: [{cw}, minE, maxE, ciLevel]\n      columns:\n"
    for c, ti, ty, ex in [('e', 'Event time', 'integer', ''), ('est', 'Estimate', 'number', ''), ('se', 'SE', 'number', ''),
                          ('lower', 'Lower', 'number', '          superTitle: Confidence interval\n'),
                          ('upper', 'Upper', 'number', '          superTitle: Confidence interval\n'),
                          ('p', 'p', 'number', '          format: zto,pvalue\n')]:
        s += f"        - name: {c}\n          title: {ti}\n          type: {ty}\n{ex}"
    s += f"\n    - name: plot\n      title: Event study\n      type: Image\n      width: 640\n      height: 420\n      renderFun: .plot\n      clearWith: [{cw}, minE, maxE, ciLevel]\n\n"
    if sp['cs'] or sp['et']:
        gcw = cw
        for nm, ti, vis, col in [('groupTab', 'Effects by cohort', 'aggGroup', ('g', 'Cohort (first treated period)')),
                                 ('calTab', 'Effects by calendar time', 'aggCalendar', ('t', 'Period'))]:
            s += (f"    - name: {nm}\n      title: {ti}\n      type: Table\n      visible: ({vis})\n      rows: 0\n      clearWith: [{gcw}, ciLevel]\n      columns:\n"
                  f"        - name: {col[0]}\n          title: {col[1]}\n          type: text\n"
                  "        - name: est\n          title: Estimate\n          type: number\n"
                  "        - name: se\n          title: SE\n          type: number\n"
                  "        - name: lower\n          title: Lower\n          type: number\n          superTitle: Confidence interval\n"
                  "        - name: upper\n          title: Upper\n          type: number\n          superTitle: Confidence interval\n\n")
        if sp['cs']: s += f"    - name: gtPlot\n      title: Group-time effects\n      type: Image\n      visible: (showGT)\n      width: 640\n      height: 420\n      renderFun: .plotGT\n      clearWith: [{gcw}]\n\n"
    s += "    - name: code\n      title: Reproducible R code\n      type: Preformatted\n      visible: (showCode)\n\n"
    s += "    - name: notes\n      title: About this estimator\n      type: Html\n      visible: (showHelp)\n      refs:\n        - roth2023\n\n"
    s += "...\n"
    return s

def b_r(sp):
    n, k = sp['name'], sp['key']
    cs = sp['cs']
    return f'''
{n}Class <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "{n}Class",
    inherit = {n}Base,
    private = list(
        .run = function() dc_single_run(self, "{k}"),
        .plot = function(image, ggtheme, theme, ...) dc_single_plot(image, ggtheme, self$options)''' + (''',
        .plotGT = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(dc_plot_gt(image$state, ggtheme, self$options$plotPalette))
            TRUE
        }''' if cs else '') + '''
    )
)
'''

for sp in SPECS:
    n = sp['name']
    open(f'jamovi/{n}.a.yaml', 'w', encoding='utf-8', newline='\n').write(a_yaml(sp))
    open(f'jamovi/{n}.u.yaml', 'w', encoding='utf-8', newline='\n').write(u_yaml(sp))
    open(f'jamovi/{n}.r.yaml', 'w', encoding='utf-8', newline='\n').write(r_yaml(sp))
    open(f'R/{n}.b.R', 'w', encoding='utf-8', newline='\n').write(b_r(sp))
print('generated', [s['name'] for s in SPECS])
