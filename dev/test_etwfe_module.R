# The analysis "Extended TWFE (Wooldridge)" run headless in the installed module (jamovi's formula sandbox is active there): linear, Poisson and logit
# models, both control groups, covariates, the cohort and calendar tables, the R code, and the messages. Run from the module root in a jamovi-like environment
# with the module installed:
#   R_LIBS_USER=notthere R_LIBS="<module R>;<jamovi>/modules/jmv/R;<jamovi>/modules/base/R" Rscript dev/test_etwfe_module.R
suppressMessages(library(didfamily))
ok <- 0L; bad <- 0L
chk <- function(what, cond) { if (isTRUE(cond)) ok <<- ok + 1L else { bad <<- bad + 1L; cat("FAIL:", what, "\n") } }
st <- read.csv("dev/csv/sim_staggered.csv"); ct <- read.csv("dev/csv/sim_etwfe.csv"); mp <- read.csv("dev/csv/mpdta.csv")

a <- didetwfe(data = st, dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", aggGroup = TRUE, aggCalendar = TRUE, showCode = TRUE, showGuide = TRUE, showHelp = TRUE)
s <- a$summary$asDF
chk("linear, not yet treated: overall ATT as in the guide (1.607)", abs(s$att - 1.607) < 0.002 && abs(s$se - 0.024) < 0.002)
bj <- didbjs(data = st, dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort")$summary$asDF
chk("linear without covariates equals the imputation estimate", abs(s$att - bj$att) < 1e-6)
chk("dynamic table with the reference row", any(a$dynamic$asDF$e == -1 & a$dynamic$asDF$est == 0) && nrow(a$dynamic$asDF) >= 4)
chk("cohort and calendar tables are filled", nrow(a$groupTab$asDF) == 3 && nrow(a$calTab$asDF) >= 5)
chk("code uses etwfe and emfx", grepl("etwfe\\(", a$code$content) && grepl("emfx", a$code$content))
chk("guide and notes present", grepl("Wooldridge", a$guide$content) && grepl("Sample data 4", a$guide$content))

n <- didetwfe(data = st, dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", etControl = "nevertreated")
chk("never treated: overall ATT 1.637 and leads in the dynamic table", abs(n$summary$asDF$att - 1.637) < 0.002 && any(n$dynamic$asDF$e == -3))

p <- didetwfe(data = ct, dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", covs = "x1", etFamily = "poisson", minE = -3, maxE = 4, aggGroup = TRUE)
chk("Poisson with a covariate: overall effect 2.60 (s.e. 0.144)", abs(p$summary$asDF$att - 2.601) < 0.003 && abs(p$summary$asDF$se - 0.144) < 0.003)
chk("Poisson: effects by cohort 2.86, 2.54, 2.14", all(abs(p$groupTab$asDF$est - c(2.862, 2.539, 2.140)) < 0.003))
chk("the method name mentions the model", grepl("poisson", p$summary$asDF$method))
l <- didetwfe(data = transform(ct, yb = as.integer(y > 4)), dep = "yb", id = "id", time = "year", treat = "g", treatType = "cohort", etFamily = "logit", etControl = "nevertreated")
chk("logit runs on a 0/1 outcome", is.finite(l$summary$asDF$att) && l$summary$asDF$att > 0)
m <- didetwfe(data = mp, dep = "lemp", id = "countyreal", time = "year", treat = "first.treat", treatType = "cohort", covs = "lpop", minE = -3, maxE = 3)
chk("mpdta with covariate lpop: overall ATT -0.051", abs(m$summary$asDF$att - (-0.0506)) < 0.002)
cl <- didetwfe(data = transform(ct, grp = id %% 50), dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", cluster = "grp", etFamily = "poisson")
chk("clustering above the unit level changes the standard error", cl$summary$asDF$se != didetwfe(data = ct, dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", etFamily = "poisson")$summary$asDF$se)

e1 <- didetwfe(data = st, dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", etFamily = "logit")$summary$asDF
chk("logit on a continuous outcome gives a message in the table", nzchar(e1$note))
e2 <- didetwfe(data = st[st$g != 0, ], dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", etControl = "nevertreated")$summary$asDF
chk("never-treated control without never-treated units gives a message", nzchar(e2$note))
e3 <- didetwfe(data = transform(ct, y = y - 6), dep = "y", id = "id", time = "year", treat = "g", treatType = "cohort", etFamily = "poisson")$summary$asDF
chk("Poisson on a negative outcome gives a message", nzchar(e3$note))
chk("the analysis has no weights option", !("weights" %in% names(formals(didetwfe))))
cat(sprintf("extended TWFE module tests: %d passed, %d failed\n", ok, bad))
if (bad) quit(status = 1)
