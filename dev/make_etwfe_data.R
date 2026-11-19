# Sample data of the extended TWFE analysis: a count outcome with staggered adoption and a covariate. Run from the module root:  Rscript dev/make_etwfe_data.R
# then Rscript dev/make_omv.R to build data/sim_etwfe.omv
set.seed(20261120)
N <- 600; Tn <- 8
d <- expand.grid(id = 1:N, year = 1:Tn)
first <- sample(c(4, 5, 6, 0), N, TRUE, c(.25, .25, .2, .3)); x1 <- round(rnorm(N), 3)
d$g <- first[d$id]; d$x1 <- x1[d$id]
treated <- d$g > 0 & d$year >= d$g
unit <- rnorm(N, 0, 0.3)[d$id]
lin <- 0.8 + 0.25 * d$x1 + 0.05 * d$year + unit + 0.15 * (d$g > 0)
eff <- treated * (0.35 + 0.08 * (d$year - d$g) + 0.15 * d$x1)
d$y <- rpois(nrow(d), exp(lin + eff))
write.csv(d[, c("id", "year", "y", "g", "x1")], "dev/csv/sim_etwfe.csv", row.names = FALSE)
cat("rows", nrow(d), " mean count", round(mean(d$y), 2), " share treated observations", round(mean(treated), 3), "\n")
