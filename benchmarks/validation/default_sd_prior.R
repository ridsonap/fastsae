# Default prior on the random-effect SD of Gaussian HB models (.default_prior_prec):
# frequentist coverage of the 95% intervals (Fay-Herriot, known vardir) for PC priors
# scaled to the response, PC(U = sd(y), alpha), and half-normal(sd(y)); EBLUP for
# reference. 15/30/60 areas (1/8 without a direct estimate), area SD 0.4/0.8,
# vardir U(0.3, 1), 150 replicates.
#
# Run from the package root: Rscript benchmarks/validation/default_sd_prior.R
#
# Result (2026-10-18), coverage / RMSE are about equal across priors:
#   m   sigma_u   EBLUP   HN(sd y)   PC a=0.01   PC a=0.05   PC a=0.10
#   15  0.4       0.960   0.966      0.917       0.942       0.950
#   30  0.4       0.940   0.958      0.905       0.929       0.937
#   60  0.4       0.927   0.943      0.898       0.918       0.926
#   15  0.8       0.930   0.932      0.853       0.890       0.901
#   30  0.8       0.934   0.929      0.888       0.909       0.917
#   60  0.8       0.942   0.939      0.919       0.928       0.931
# The PC priors put much mass near sigma = 0, shrink the area SD and give too short
# intervals; the half-normal prior covers as the EBLUP intervals. The unscaled PC(1, 0.01)
# prior (INLA's default) also depends on the unit of y: covering 94%, 85% and 38% with
# y multiplied by 0.1, 1 and 10 (30 areas, 60 replicates).
devtools::load_all(quiet = TRUE)
q <- function(e) suppressWarnings(suppressMessages(e))
R <- 150
out <- NULL
set.seed(10)
for (m in c(15, 30, 60)) {
  for (su in c(0.4, 0.8)) {
    x <- rnorm(m)
    psi <- runif(m, 0.3, 1)
    miss <- seq_len(max(2, round(m / 8)))
    for (r in seq_len(R)) {
      theta <- 1 + x + rnorm(m, sd = su)
      y <- theta + rnorm(m, sd = sqrt(psi))
      y[miss] <- NA
      d <- data.frame(y = y, x = x, v = psi)
      sy <- sd(y, na.rm = TRUE)
      pri <- list(
        "PC a=0.01" = list(prior = "pc.prec", param = c(sy, 0.01)),
        "PC a=0.05" = list(prior = "pc.prec", param = c(sy, 0.05)),
        "PC a=0.10" = list(prior = "pc.prec", param = c(sy, 0.10)),
        "HN(sd y)" = .half_normal_prec(sy)
      )
      for (nm in names(pri)) {
        f <- tryCatch(q(hb_area(y ~ x,
          data = d, vardir = "v", print_result = FALSE,
          engine = inla_engine(prior_prec = pri[[nm]])
        )), error = function(e) NULL)
        if (is.null(f)) next
        h <- f$df_hb
        out <- rbind(out, data.frame(
          m = m, sigma_u = su, prior = nm, unsamp = seq_len(m) %in% miss,
          cover = theta >= h$ci_lower & theta <= h$ci_upper, w = h$ci_upper - h$ci_lower,
          err2 = (h$hb - theta)^2
        ))
      }
      fe <- q(eblup_fh(y ~ x, vardir = ~v, data = d, print_result = FALSE))$df_eblup
      out <- rbind(out, data.frame(
        m = m, sigma_u = su, prior = "EBLUP", unsamp = seq_len(m) %in% miss,
        cover = abs(fe$eblup - theta) <= 1.96 * sqrt(fe$mse), w = 2 * 1.96 * sqrt(fe$mse), err2 = (fe$eblup - theta)^2
      ))
    }
  }
}
a <- aggregate(cbind(cover, w, err2) ~ m + sigma_u + prior, out, mean)
a$rmse <- sqrt(a$err2)
a$err2 <- NULL
print(a[order(a$sigma_u, a$m, a$prior), ], digits = 3, row.names = FALSE)
