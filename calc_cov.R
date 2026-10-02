fl <- getOption("covr.flags")
if (is.null(fl)) fl <- character()
fl["FLIBS"] <- "-L/opt/homebrew/Cellar/gcc/16.2.0/lib/gcc/16 -L/opt/homebrew/Cellar/gcc/16.2.0/lib/gcc/current/gcc/aarch64-apple-darwin25/16 -lgcov"
options(covr.flags = fl)

cat("Starting covr::package_coverage()...\n")
flush.console()

cv <- covr::package_coverage(quiet = FALSE)
saveRDS(cv, "cov_res.rds")
print(cv)
pct <- covr::percent_coverage(cv)
cat(sprintf("\n=== OVERALL COVERAGE: %.2f%% ===\n", pct))
