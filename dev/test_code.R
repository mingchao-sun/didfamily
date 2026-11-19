# Checks that the R code generated for users reproduces the module's own results.
# Runs in PLAIN R (no jamovi sandbox): the generated script uses sunab() and i().
#   R_LIBS_USER=notthere R_LIBS="<module>/R;<jamovi>/Resources/modules/jmv/R;<jamovi>/Resources/modules/base/R" Rscript dev/test_code.R
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_estimators", "core_did2s", "vendor_honestdid", "core_code")) source(file.path(root, "R", paste0(f, ".R")))
suppressMessages({ library(fixest); library(did); library(didimputation); library(did2s) })

check <- function(cond, msg) {
    cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
    if (!isTRUE(cond)) quit(status = 1)
}

run_case <- function(label, dat, y, id, time, treat, type, covs = NULL, cluster = NULL, weights = NULL, informational = character()) {
    cat("\n---", label, "---\n")
    d <- dc_prepare(dat, y, id, time, treat, type, covs = covs, cluster = cluster, weights = weights)
    o <- list(dep = y, id = id, time = time, treat = treat, treatType = type, covs = covs, cluster = cluster, weights = weights,
              minE = -4, maxE = 4, twfe = TRUE, cs = TRUE, sa = TRUE, bjs = TRUE, gardner = TRUE, stacked = TRUE,
              csControl = "nevertreated", csMethod = "dr", boot = FALSE, bootReps = 999, seed = 1234)
    code <- dc_code_compare(o, d$info)
    code <- sub("read.csv(\"your_data.csv\")", "dat_test", code, fixed = TRUE)
    env <- new.env(parent = globalenv()); env$dat_test <- dat
    out <- tryCatch({ suppressMessages(eval(parse(text = code), env)); NULL }, error = function(e) conditionMessage(e))
    check(is.null(out), paste("generated script runs without error", if (!is.null(out)) out else ""))
    tab <- env$summary_table
    mod <- list(TWFE = dc_twfe(d, -4, 4), `Callaway-Sant'Anna` = dc_cs(d, -4, 4, "nevertreated", "dr"),
                `Sun-Abraham` = dc_sa(d, -4, 4), `Borusyak-Jaravel-Spiess` = dc_bjs(d, -4, 4),
                `Gardner (did2s)` = dc_gardner(d, -4, 4), `Stacked DiD` = dc_stacked(d, -4, 4))
    for (nm in names(mod)) {
        if (!isTRUE(mod[[nm]]$ok)) {
            check(nrow(tab[tab$estimator == nm, ]) == 0, paste(nm, ": skipped in the module and in the script"))
            next
        }
        r <- tab[tab$estimator == nm, ]
        if (nm %in% informational) {
            cat(sprintf("  %-26s module %.6f  script %.6f  (informational: aggregation weights differ by design)\n",
                        nm, mod[[nm]]$att, r$att))
            next
        }
        cat(sprintf("  %-26s module %.6f (%.5f)  script %.6f (%.5f)\n", nm, mod[[nm]]$att, mod[[nm]]$se, r$att, r$se))
        check(abs(r$att - mod[[nm]]$att) < 1e-6, paste(nm, ": script ATT equals module ATT"))
        check(abs(r$se - mod[[nm]]$se) < 1e-5, paste(nm, ": script SE equals module SE"))
    }
    invisible(list(env = env, mod = mod, code = code))
}

sim <- read.csv(file.path(root, "dev/csv/sim_staggered.csv"))
res1 <- run_case("indicator, no covariates", sim, "y", "id", "year", "D", "indicator")
# event-study pieces
e <- res1$env; m <- res1$mod
cs_es <- e$cs_es; i <- which(cs_es$egt == 2)
check(abs(cs_es$att.egt[i] - m$`Callaway-Sant'Anna`$dynamic$est[m$`Callaway-Sant'Anna`$dynamic$e == 2]) < 1e-8, "script CS event-time estimate equals module (e = 2)")
tw_es <- coeftable(e$twfe_es)
check(abs(tw_es["dc_rel::3", 1] - m$TWFE$dynamic$est[m$TWFE$dynamic$e == 3]) < 1e-8, "script TWFE event-time estimate equals module (e = 3)")
bes <- as.data.frame(e$bjs_es); bm <- m$`Borusyak-Jaravel-Spiess`$dynamic
for (k in c(-3, 0, 2)) {
    sc <- bes$estimate[bes$term == as.character(k)]
    check(length(sc) == 1 && abs(sc - bm$est[bm$e == k]) < 1e-6, sprintf("script BJS event-time estimate equals module (e = %d)", k))
}

# cohort-type input, covariates (numeric + factor), cluster variable
set.seed(5); N2 <- 600; T2 <- 8
x <- rnorm(N2); grp <- sample(c("a", "b", "c"), N2, TRUE); cl <- sample(1:60, N2, TRUE)
g2 <- ifelse(runif(N2) < .65, ifelse(runif(N2) < plogis(1.5 * x), 4, 6), 0)
cv <- expand.grid(id = 1:N2, year = 1:T2)
cv$g <- g2[cv$id]; cv$x <- x[cv$id]; cv$grp <- grp[cv$id]; cv$cl <- cl[cv$id]; cv$z <- rnorm(nrow(cv))
cv$D <- as.integer(cv$g > 0 & cv$year >= cv$g)
cv$y <- rnorm(N2)[cv$id] + rnorm(T2, sd = .3)[cv$year] + 0.4 * cv$x * cv$year + 0.5 * cv$z + cv$D + rnorm(nrow(cv), sd = .5)
run_case("first-treated-period input, covariates x, z, grp (factor)", cv, "y", "id", "year", "g", "cohort", covs = c("x", "z", "grp"))

# cluster variable (forces the CS bootstrap) + equally spaced periods: the script must still reproduce the module
cv2 <- cv
cv2$year <- 2000 + 2 * (cv2$year - 1)
cv2$g <- ifelse(cv2$g > 0, 2000 + 2 * (cv2$g - 1), 0)
run_case("cluster variable + bootstrap + equally spaced periods", cv2, "y", "id", "year", "g", "cohort", covs = "x", cluster = "cl")
# negative periods, indicator input
neg <- sim; neg$year <- neg$year - 6
run_case("negative periods, indicator", neg, "y", "id", "year", "D", "indicator")
# sampling weights: the script must reproduce TWFE, CS, Gardner and stacked; BJS is skipped in both
wcv <- cv
wcv$pop <- 1 + (wcv$id %% 5)
run_case("sampling weights", wcv, "y", "id", "year", "g", "cohort", covs = "x", weights = "pop")
cat("\nAll generated-code tests passed.\n")
