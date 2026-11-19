# M0 feasibility: can DiD engine packages be installed/loaded in jamovi's bundled R?
# Run with: "<jamovi>/Frameworks/R/bin/Rscript.exe" dev/check_deps.R
pkgs <- c("fixest","did","DRDID","DIDmultiplegtDYN","didimputation","did2s",
          "etwfe","bacondecomp","TwoWayFEWeights","fect","triplediff",
          "HonestDiD","pretrends","synthdid","contdid")
repos <- "https://cloud.r-project.org"
res <- data.frame(pkg = pkgs, installed = NA, loads = NA, version = NA, note = "")
for (i in seq_along(pkgs)) {
  p <- pkgs[i]
  if (!requireNamespace(p, quietly = TRUE)) {
    av <- tryCatch(rownames(available.packages(repos = repos)), error = function(e) character())
    if (p %in% av) {
      tryCatch(install.packages(p, repos = repos, quiet = TRUE),
               error = function(e) res$note[i] <<- conditionMessage(e))
    } else res$note[i] <- "not on CRAN (GitHub-only?)"
  }
  res$installed[i] <- requireNamespace(p, quietly = TRUE)
  res$loads[i] <- tryCatch({ suppressMessages(loadNamespace(p)); TRUE }, error = function(e) FALSE)
  if (res$installed[i]) res$version[i] <- as.character(packageVersion(p))
}
print(res, row.names = FALSE)
