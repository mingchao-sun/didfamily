# Classic-data check: Castle Doctrine laws (Cheng & Hoekstra), shipped with bacondecomp.
# Staggered adoption 2005-2009, 29 never-treated states, balanced 2000-2010 panel.
#   R_LIBS_USER=notthere R_LIBS="<module>/R;<jamovi>/Resources/modules/jmv/R;<jamovi>/Resources/modules/base/R" Rscript dev/test_castle.R
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_estimators", "core_did2s", "vendor_honestdid", "core_diag")) source(file.path(root, "R", paste0(f, ".R")))
source(file.path(root, "dev/sandbox.R")); dc_enable_sandbox()
check <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) quit(status = 1) }

castle <- as.data.frame(bacondecomp::castle)
cat("columns:", paste(names(castle)[1:8], collapse = ", "), "...\n")
d <- dc_prepare(castle, "l_homicide", "state", "year", "post", "indicator")
i <- d$info
cat(sprintf("%d states, %d years, %d never treated, cohorts: %s\n", i$n_units, i$n_periods, i$n_never,
            paste(sprintf("%s (n=%d)", names(i$cohorts), as.integer(i$cohorts)), collapse = ", ")))
check(i$n_units == 50 && i$n_periods == 11 && i$balanced, "50 states, 11 years, balanced")
check(i$n_never == 29 && sum(i$cohorts) == 21, "29 never-treated and 21 treated states")

# TWFE: the headline Cheng-Hoekstra estimate is about +8-10% homicide
tw <- dc_twfe(d, -4, 4)
cat(sprintf("TWFE ATT %.4f (se %.4f)\n", tw$att, tw$se))
check(tw$ok && tw$att > 0.05 && tw$att < 0.12, "TWFE estimate is positive and close to the published ~0.08-0.10")

# Goodman-Bacon decomposition reproduces the TWFE coefficient exactly
bc <- dc_bacon(d)
print(bc$summary)
check(bc$ok && abs(bc$twfe - tw$att) < 1e-6, "Bacon decomposition reproduces the TWFE coefficient")
check(abs(sum(bc$detail$weight) - 1) < 1e-8, "Bacon weights sum to one")
x <- bacondecomp::bacon(l_homicide ~ post, data = castle, id_var = "state", time_var = "year", quietly = TRUE)
check(abs(sum(x$weight * x$estimate) - bc$twfe) < 1e-10, "identical to a direct bacondecomp::bacon call")
w <- dc_twfe_weights(d); cat(sprintf("negative TWFE weights on %d of %d treated obs\n", w$neg_n, w$n_treated))
check(w$ok && abs(sum(w$w) - 1) < 1e-8, "TWFE weights sum to one")

# the six estimators agree on sign and rough size; CS and imputation need not equal TWFE
res <- list(tw, dc_cs(d, -4, 4), dc_sa(d, -4, 4), dc_bjs(d, -4, 4), dc_gardner(d, -4, 4), dc_stacked(d, -4, 4))
for (r in res) cat(sprintf("  %-26s ATT %.4f (se %.4f) %s\n", r$method, r$att, r$se, r$message))
check(all(vapply(res, function(r) isTRUE(r$ok), TRUE)), "all six estimators run on the Castle data")
check(all(vapply(res, function(r) r$att > 0, TRUE)), "all six estimators find a positive effect on homicide")
bjs <- res[[4]]; gd <- res[[5]]
check(abs(bjs$att - gd$att) < 1e-6, "BJS and Gardner agree on the Castle data (both imputation)")
# pre-trend tests run
pt <- dc_pretests(d, -4)
print(pt[, 1:4])
check(all(is.finite(pt$p[pt$method != "Callaway-Sant'Anna"])), "TWFE and Sun-Abraham pre-trend tests return p-values on the Castle data")
cs_row <- pt[pt$method == "Callaway-Sant'Anna", ]
check(is.finite(cs_row$p) || nzchar(cs_row$note), "CS pre-test either returns a p-value or explains why not (small cohorts make did's Wald matrix singular)")
cat("\nAll Castle-data checks passed.\n")

# population-weighted analysis (weights = state population in 2000), as in the applied literature
cat("\n--- weighted Castle ---\n")
castle_w <- read.csv(file.path(root, "dev/csv/castle.csv"))
dw <- dc_prepare(castle_w, "l_homicide", "state", "year", "post", "indicator", weights = "pop2000")
rw <- list(dc_twfe(dw, -4, 4), dc_cs(dw, -4, 4), dc_sa(dw, -4, 4), dc_gardner(dw, -4, 4), dc_stacked(dw, -4, 4))
for (r in rw) cat(sprintf("  %-26s ATT %.4f (se %.4f, N=%d) %s\n", r$method, r$att, r$se, r$n_obs, r$message))
check(all(vapply(rw, function(r) isTRUE(r$ok), TRUE)), "all weighted estimators run on the Castle data")
check(abs(rw[[1]]$att - tw$att) > 1e-4, "population weights change the TWFE estimate")
check(rw[[1]]$df == 49, "TWFE uses t(49) with 50 states")
check(!dc_bjs(dw, -4, 4)$ok, "BJS is skipped with weights")
cat("\nAll weighted Castle checks passed.\n")
