# Tests for the classic 2x2 DiD core (run under the emulated jamovi sandbox).
root <- if (file.exists("R/core_prepare.R")) "." else ".."
for (f in c("core_prepare", "core_estimators", "core_2x2")) source(file.path(root, "R", paste0(f, ".R")))
source(file.path(root, "dev/sandbox.R")); dc_enable_sandbox()
check <- function(cond, msg) {
    cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
    if (!isTRUE(cond)) quit(status = 1)
}

set.seed(2); N <- 2000
x <- rnorm(N); trt <- rbinom(N, 1, plogis(0.8 * x)); grp <- sample(c("a", "b"), N, TRUE)
long <- data.frame(id = rep(1:N, 2), year = rep(c(2010, 2015), each = N), trt = rep(trt, 2), x = rep(x, 2), grp = rep(grp, 2))
long$y <- 1 + 0.5 * long$x * (long$year == 2015) + 0.3 * (long$year == 2015) + 0.2 * long$trt +
          1.5 * long$trt * (long$year == 2015) + rnorm(2 * N, sd = .5)       # true effect 1.5; trends depend on x

d0 <- dc_prepare22(long, "y", "trt", "year", id = "id")
check(d0$panel && d0$labels[2] == "2015" && d0$info$n_units == N, "prepare: panel, post period = later value")
r0 <- dc_did22(d0, c("ols"))
check(r0$ols$ok && abs(r0$ols$att - 1.5) < 0.5, "no covariates: OLS DiD runs")

dx <- dc_prepare22(long, "y", "trt", "year", id = "id", covs = "x")
rx <- dc_did22(dx)
for (m in names(rx)) cat(sprintf("  %-30s %.3f (se %.3f) %s\n", rx[[m]]$method, rx[[m]]$att, rx[[m]]$se, rx[[m]]$message))
check(all(vapply(rx, function(r) r$ok, TRUE)), "all five estimators run under the sandbox")
check(abs(rx$dr$att - 1.5) < 0.12 && abs(rx$drtrad$att - 1.5) < 0.12, "doubly robust estimators recover the effect")
check(abs(rx$ols$att - 1.5) > 0.2, "plain OLS with a linear covariate stays biased (trends depend on x x post)")
check(abs(rx$ipw$att - 1.5) < 0.2 && abs(rx$or$att - 1.5) < 0.2, "IPW and outcome regression are close to the truth")

# repeated cross-sections give nearly the same answer
dr <- dc_prepare22(long, "y", "trt", "year", covs = "x")
check(!dr$panel, "no unit ID -> repeated cross-sections")
rr <- dc_did22(dr)
check(all(vapply(rr, function(r) r$ok, TRUE)) && abs(rr$dr$att - rx$dr$att) < 0.05, "cross-section DR estimate close to panel DR estimate")

# factor covariate, period given as text
lf <- long; lf$period <- ifelse(lf$year == 2015, "after", "before")
dfc <- dc_prepare22(lf, "y", "trt", "period", id = "id", covs = c("x", "grp"))
check(dfc$labels[2] == "before" || dfc$labels[2] == "after", "text period variable accepted")
check(nrow(dfc$covs) == 2, "factor covariate becomes a dummy")
check(all(vapply(dc_did22(dfc), function(r) r$ok, TRUE)), "all estimators run with numeric and factor covariates")

m4 <- dc_means22(d0)
check(nrow(m4) == 4 && all(m4$n == N / 1 * c(1, 1, 1, 1) * 0 + m4$n) && sum(m4$n) == 2 * N, "cell means table has four cells covering all rows")
did_by_hand <- (m4$mean[m4$group == 1 & m4$post == 1] - m4$mean[m4$group == 1 & m4$post == 0]) -
               (m4$mean[m4$group == 0 & m4$post == 1] - m4$mean[m4$group == 0 & m4$post == 0])
check(abs(did_by_hand - r0$ols$att) < 1e-8, "OLS interaction equals the difference of difference in means")

# errors
bad3 <- long; bad3$year[1:5] <- 2020
check(inherits(try(dc_prepare22(bad3, "y", "trt", "year", id = "id"), silent = TRUE), "try-error"), "three periods rejected")
bad_g <- long; bad_g$trt[1] <- 2
check(inherits(try(dc_prepare22(bad_g, "y", "trt", "year", id = "id"), silent = TRUE), "try-error"), "group values other than 0/1 rejected")
unb <- long[-1, ]
check(inherits(try(dc_prepare22(unb, "y", "trt", "year", id = "id"), silent = TRUE), "try-error"), "unbalanced two-period panel rejected with a unit ID")
check(!inherits(try(dc_prepare22(unb, "y", "trt", "year"), silent = TRUE), "try-error"), "the same data are fine as repeated cross-sections")
cat("\nAll 2x2 tests passed.\n")
tg <- long; tg$trt <- ifelse(tg$trt == 1, "treated", "control")
check(grepl("numeric 0/1", tryCatch(dc_prepare22(tg, "y", "trt", "year", id = "id"), error = function(e) conditionMessage(e))), "text group variable: clear message")
