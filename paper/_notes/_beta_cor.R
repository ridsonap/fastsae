x <- readRDS("/Volumes/work/_MainR/fastsae-cov/inst/extdata/beta_estimates_n1000.rds")
a <- x[["estimasi fastsae"]]; b <- x[["estimasi tipsae"]]
ma <- x[["mse fastsae"]]; mb <- x[["mse tipsae"]]
cat("pearson estimasi:", sprintf("%.7f", cor(a, b)), "\n")
cat("spearman estimasi:", sprintf("%.7f", cor(a, b, method = "spearman")), "\n")
cat("MAE:", sprintf("%.7f", mean(abs(a - b))), "\n")
cat("RMSD:", sprintf("%.7f", sqrt(mean((a - b)^2))), "\n")
cat("mse pearson:", sprintf("%.7f", cor(ma, mb)), "\n")
cat("mse MAD:", sprintf("%.7f", mean(abs(ma - mb))), "\n")
