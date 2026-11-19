# Timing smoke test on a large panel (jamovi-like library path, sandbox emulated).
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_estimators", "core_diag", "core_code")) source(file.path(root, "R", paste0(f, ".R")))
source(file.path(root, "dev/sandbox.R")); dc_enable_sandbox()
N <- as.integer(commandArgs(TRUE)[1]); T <- as.integer(commandArgs(TRUE)[2])
set.seed(1); nc <- 8
cohorts <- sort(sample(4:(T - 2), nc)); g <- sample(c(cohorts, 0), N, TRUE)
d <- expand.grid(id = 1:N, year = 1:T); d$g <- g[d$id]; d$D <- as.integer(d$g > 0 & d$year >= d$g)
d$y <- rnorm(N)[d$id] + 0.1 * d$year + d$D * (0.5 + 0.1 * (d$year - d$g)) + rnorm(nrow(d))
cat(sprintf("panel: %d units x %d periods = %d rows, %d cohorts\n", N, T, nrow(d), nc))
tm <- function(label, expr) { t0 <- Sys.time(); r <- force(expr); cat(sprintf("  %-22s %6.1f s  %s\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs")), if (is.list(r) && !is.null(r$ok)) if (r$ok) "ok" else r$message else "")); invisible(r) }
which <- commandArgs(TRUE)[3]
p <- tm("prepare", dc_prepare(d, "y", "id", "year", "D", "indicator"))
fns <- list(twfe = function() dc_twfe(p, -4, 4), cs = function() dc_cs(p, -4, 4), sa = function() dc_sa(p, -4, 4),
            bjs = function() dc_bjs(p, -4, 4), gardner = function() dc_gardner(p, -4, 4), stacked = function() dc_stacked(p, -4, 4),
            bacon = function() dc_bacon(p), weights = function() dc_twfe_weights(p), pretests = function() dc_pretests(p, -4))
if (is.na(which)) which <- names(fns)
for (w in which) { cat("running", w, "...
"); flush.console(); tm(w, fns[[w]]()); flush.console() }
