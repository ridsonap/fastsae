library(devtools)


# chmod +x configure
Rcpp::compileAttributes()
devtools::document()
devtools::check()
devtools::test()
load_all()




# SFH ---------------------------------------------------------------------
# Spatial Fay-Herriot model
m1 <- eblup_sfh(
  y ~ x1 + x2 + x3,
  data = mys,
  domain = ~area,
  vardir = ~vardir,
  W = mys_proxmat
)

# Spatial Fay-Herriot model with Parametric Bootstrap MSE
m2 <- eblup_sfh(
  y ~ x1 + x2 + x3,
  data = mys,
  domain = ~area,
  vardir = ~vardir,
  mse_method = "pbmse",
  B = 50,
  W = mys_proxmat
)


# BHF ---------------------------------------------------------------------
library(dplyr)
df_meanpop <- cornsoybeanmeans |>
  rename(CornPix = MeanCornPixPerSeg, SoyBeansPix = MeanSoyBeansPixPerSeg)
df_cornsoybean <- cornsoybean |>
  rename(CountyIndex = County)

res <- eblup_bhf(
  formula = CornHec ~ CornPix + SoyBeansPix,
  Xpop = df_meanpop,
  unit_data = df_cornsoybean,
  domain_var = "CountyIndex",
  popsize_var = "PopnSegments"
)



# EBP Area-Level (INLA & Laplace) -----------------------------------------
# 1. Non-spatial Fay-Herriot via INLA
m_ebp0 <- ebp_area(
  y ~ x1 + x2 + x3,
  data = mys,
  domain = ~area,
  vardir = ~vardir,
  family = "gaussian"
)

# 2. Spatial BYM2 Fay-Herriot via INLA
m_ebp_bym2 <- ebp_area(
  y ~ x1 + x2 + x3,
  data = mys,
  domain = ~area,
  vardir = ~vardir,
  W = mys_proxmat,
  spatial = "bym2"
)

# 3. Model comparison visualization
autoplot(list("FH-INLA (Non-spatial)" = m_ebp0, "SFH-INLA (BYM2)" = m_ebp_bym2), type = "comparison")

# 4. Diagnostic evaluation of estimates
diag_fh <- diagnose(m_ebp0)
print(diag_fh)

diag_sfh <- diagnose(m_ebp_bym2)
print(diag_sfh)
autoplot(diag_sfh, type = "calibration")


# Simulation & Multi-Distribution Area Data -------------------------------
# 1. Generate spatial proximity matrix (KNN, Grid, Ring)
W_knn <- sim_spatial_weights(D = 40, type = "knn", k = 4, style = "B", seed = 123)
W_grid <- sim_spatial_weights(D = 36, type = "grid", style = "W")

# 2. Simulate area-level data with multi-distribution responses & spatial effects
sim_res <- sim_area_data(D = 42, spatial_type = "knn", rho = 0.5, phi = 0.6, n_unsampled = 6, seed = 2026)
print(sim_res)
head(sim_res$data)

# 3. Fit Poisson EBP with INLA BYM2 spatial model
fit_pois <- ebp_area(
  y_poisson ~ x1 + x2,
  data = sim_res$data,
  exposure = "exposure",
  family = "poisson",
  W = sim_res$W,
  spatial = "bym2"
)
summary(fit_pois)

# 4. Using built-in sim_area with mys_proxmat
data(sim_area)
data(mys_proxmat)
fit_multi_bin <- ebp_area(
  y_binomial ~ x1 + x2,
  data = sim_area,
  trials = "trials",
  family = "binomial",
  W = mys_proxmat,
  spatial = "bym2"
)
summary(fit_multi_bin)


# Spatio-Temporal Series Simulation ---------------------------------------
# 1. Simulate 30 domains over 4 years with SAR spatial and AR(1) temporal dynamics
sim_panel_gen <- sim_series_data(
  D = 30,
  T = 4,
  time_start = 2021,
  rho_s = 0.5,
  rho_t = 0.6,
  n_unsampled = 0,
  prop_intermittent = 0,
  seed = 123
)
print(sim_panel_gen)
head(sim_panel_gen$data)

# 2. Fit Spatio-Temporal Fay-Herriot Model (eblup_stfh)
fit_stfh <- eblup_stfh(
  y_gaussian ~ x1 + x2,
  data = sim_panel_gen$data,
  domain = ~area,
  time = ~year,
  vardir = ~vardir,
  W = sim_panel_gen$W_std
)
summary(fit_stfh)

# 3. Using built-in sim_panel dataset
data(sim_panel)
head(sim_panel)



# Plot -------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

eblup <- readRDS("inst/exdata/eblup_benchmark.rds")
seblup <- readRDS("inst/exdata/seblup_benchmark.rds")
steblup <- readRDS("inst/exdata/steblup_benchmark.rds")

eblup_time <- eblup %>% mutate(Algorithm = "Fay Herriot")
seblup_time <- seblup %>% mutate(Algorithm = "Spatial FH")
steblup_time <- steblup %>% mutate(Algorithm = "Spatio Temporal FH")


combined_data <- bind_rows(eblup_time, seblup_time, steblup_time) %>%
  rename(Median = `Median (s)`, Mem = `Memory (MB)`)

ggplot(combined_data, aes(x = n, y = Median, color = Method)) +
  geom_line(aes(alpha = Method), linewidth = 1) +
  geom_point(aes(alpha = Method), size = 2.5) +
  scale_alpha_manual(
    NULL,
    values = c(0.3, 1, 0.3),
    guide = 'none'
  ) +
  scale_x_log10(breaks = c(30, 50, 100, 250, 500, 1000)) +
  scale_y_log10(labels = label_number(suffix = "s")) +
  facet_wrap(~Algorithm) +
  labs(
    title = "Execution Time Comparison (Log Scale)",
    subtitle = "Lower is faster",
    x = "Number of Domains (n)",
    y = "Median Time (seconds)",
    color = "Package"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  'README_files/figure-gfm/benchmark_plots-1.png',
  scale = 0.9,
  dpi = 300,
  width = 8,
  height = 4
)


ggplot(combined_data, aes(x = n, y = Mem, color = Method)) +
  geom_line(aes(alpha = Method), linewidth = 1) +
  geom_point(aes(alpha = Method), size = 2.5) +
  scale_alpha_manual(
    NULL,
    values = c(0.3, 1, 0.3),
    guide = 'none'
  ) +
  scale_x_log10(breaks = c(30, 50, 100, 250, 500, 1000)) +
  scale_y_log10(labels = label_number(suffix = "MB")) +
  facet_wrap(~Algorithm) +
  labs(
    title = "Memory Usage Comparison (Log Scale)",
    subtitle = "Lower is better",
    x = "Number of Domains (n)",
    y = "Memory usage (MB)",
    color = "Package"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  'README_files/figure-gfm/benchmark_plots_mem.png',
  scale = 0.9,
  dpi = 300,
  width = 8,
  height = 4
)


combined_data |>
  filter(Algorithm == "Spatio Temporal FH") |>
  group_by(Method) |>
  summarise(
    median = mean(Median),
    mem = mean(Mem)
  )


dim(mys_proxmat)
?sae::eblupSFH()

