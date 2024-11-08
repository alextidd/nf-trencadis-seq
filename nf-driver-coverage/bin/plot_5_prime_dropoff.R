#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(ggplot2)
library(optparse)

# parse arguments
option_list <- list(
  make_option("--gene", type = "character"),
  make_option("--min_cov", type = "numeric"),
  make_option("--plot_device", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
saveRDS(opts, "opts.rds")

# read in exonic coverage
exonic_cov <-
  list.files(pattern = "_coverage_per_exonic_position.tsv") %>%
  purrr::map(function(file) {
    run <- paste(strsplit(file, "_")[[1]][1:3], collapse = "_")
    readr::read_tsv(file) %>%
      dplyr::mutate(run = run)
  }) %>%
  dplyr::bind_rows()

p <-
  exonic_cov %>%
  dplyr::arrange(dplyr::desc(coverage)) %>%
  ggplot(aes(x = exonic_distance_from_5_prime, y = n_cells)) +
  # plot grey background (cov = 0)
  geom_col(fill = NA, colour = "grey86", position = "stack", width = 1,
           linewidth = 0.0001) +
  geom_col(fill = "grey86", position = "stack", width = 1) +
  # plot dark grey background (cov > 0) +
  geom_col(data = . %>% dplyr::filter(coverage > 0),
           fill = NA, colour = "grey48", position = "stack", width = 1,
           linewidth = 0.0001) +
  geom_col(data = . %>% dplyr::filter(coverage > 0),
           fill = "grey48", position = "stack", width = 1) +
  # plot coverage (cov >= min_cov)
  geom_col(data = . %>% dplyr::filter(coverage >= opts$min_cov),
           aes(colour = coverage), position = "stack", width = 1,
           linewidth = 0.0001) +
  geom_col(data = . %>% dplyr::filter(coverage >= opts$min_cov),
           aes(fill = coverage), position = "stack", width = 1) +
  viridis::scale_fill_viridis(na.value = "grey") +
  viridis::scale_colour_viridis(na.value = "grey") +
  facet_grid(run ~ ., scales = "free_y") +
  theme_minimal()

# save
ggsave(paste0(opts$gene, "_exonic_coverage_plot.", opts$plot_device),
       p, dpi = 500)
saveRDS(p, paste0(opts$gene, "_exonic_coverage_plot.rds"))