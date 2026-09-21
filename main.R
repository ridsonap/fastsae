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

