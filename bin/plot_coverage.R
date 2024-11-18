#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(ggplot2)
library(patchwork)
library(optparse)

# parse arguments
option_list <- list(
  make_option("--meta_run", type = "character"),
  make_option("--meta_seq_type", type = "character"),
  make_option("--meta_kit", type = "character"),
  make_option("--gene", type = "character"),
  make_option("--min_cov", type = "numeric"),
  make_option("--plot_device", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# read plot data
g <- readRDS(paste0(opts$meta_run, "_", opts$gene, ".rds"))
regions <- readRDS(paste0(opts$meta_run, "_", opts$gene, "_regions", ".rds"))
features <- readRDS(paste0(opts$meta_run, "_", opts$gene, "_features.rds"))
mut_annots <- readRDS(paste0(opts$meta_run, "_", opts$gene, "_mut_annots.rds"))
cov_per_ct <- readRDS(paste0(opts$meta_run, "_", opts$gene, "_cov_per_ct.rds"))

# functions
p_facet_theme <-
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "grey", fill = NA,
                                    linewidth = 0.1),
        strip.background = element_rect(color = "grey", fill = NA,
                                        linewidth = 0.1, linetype = "solid"))

plot_wgs_mutations <-
  function(p_dat, facet_by_region = FALSE,
           muts_subtitle = "SNVs called in WGS from Mitchell 2022") {
    n_muts <-
      p_dat %>%
      dplyr::filter(!is.na(mut_ref)) %>%
      dplyr::pull(mut_ref) %>%
      dplyr::n_distinct()
    p <-
      p_dat %>%
      dplyr::select(dplyr::any_of(c("pos", "mut_type", "ref", "alt",
                                    "region"))) %>%
      dplyr::mutate(mut = paste(ref, alt, sep = ">")) %>%
      dplyr::distinct() %>%
      ggplot(aes(x = pos, y = 0, colour = mut_type)) +
      geom_point(na.rm = TRUE) +
      theme_void() +
      scale_colour_discrete(na.translate = FALSE, na.value = NA) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      labs(subtitle = paste(n_muts, muts_subtitle))
    if (facet_by_region) {
      p <-
        p +
        facet_grid(. ~ region, scales = "free_x", space = "free_x") +
        p_facet_theme
    }
    p
  }

plot_mutations <-
  function(p_dat, facet_by_region = FALSE,
           muts_subtitle = "SNVs called in WGS from Mitchell 2022") {
    p <-
      p_dat %>%
      dplyr::mutate(
        mut_type = dplyr::case_when(alt_depth == 0 | is.na(alt_depth) ~ NA,
                                    TRUE ~ mut_type),
        mut_ref = dplyr::case_when(alt_depth == 0 | is.na(alt_depth) ~ NA,
                                   TRUE ~ mut_ref)) %>%
      plot_wgs_mutations(facet_by_region = facet_by_region,
                         muts_subtitle = muts_subtitle)
  }

plot_transcripts <- function(p_dat, facet_by_region = FALSE) {
  p <-
    p_dat %>%
    ggplot(aes(x = start, xend = end,
                y = transcript_id, yend = transcript_id, size = type)) +
    geom_segment(data = . %>% dplyr::filter(type != "start_codon"),
                  aes(colour = exon_alternating)) +
    scale_size_manual(values = c("transcript" = 0.5, "exon" = 3)) +
    scale_colour_manual(values = c("1" = "royalblue4", "0" = "steelblue3")) +
    facet_grid(transcript_type ~ ., space = "free_y", scales = "free_y") +
    theme_minimal() +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          strip.text.x = element_blank(),
          axis.text.x = element_text(angle = -90, hjust = 1),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(colour = "none") +
    labs(x = "position")
  if (facet_by_region) {
    p <-
      p +
      facet_grid(transcript_type ~ region, space = "free", scales = "free") +
      scale_x_discrete(expand = c(0, 0)) +
      p_facet_theme +
      theme(panel.spacing.y = unit(0.5, "lines"))
  }
  p
}

plot_cov <- function(p_dat, facet_by_region = FALSE, y_var = "n_cells",
                     p_title = "Coverage per cell") {
  if (y_var == "percent_cells") {
    p_dat <-
      p_dat %>%
      dplyr::filter(coverage > 0) %>%
      dplyr::mutate(percent_cells = 100 * n_cells / total_cells)
  }
  p <-
    p_dat %>%
    dplyr::arrange(-coverage) %>%
    ggplot(aes(x = pos, y = get(y_var))) +
    # plot grey background (cov = 0)
    geom_col(fill = NA, colour = "grey86", position = "stack", width = 1,
             linewidth = 0.0001) +
    geom_col(fill = "grey86", position = "stack", width = 1) +
    # plot dark grey background (cov > 0)
    geom_col(data = . %>% dplyr::filter(coverage > 0),
             fill = NA, colour = "grey48", position = "stack", width = 1,
             linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage > 0),
             fill = "grey48", position = "stack", width = 1) +
    # plot coverage
    geom_col(data = . %>% dplyr::filter(coverage >= opts$min_cov),
             aes(colour = coverage), position = "stack", width = 1,
             linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage >= opts$min_cov),
             aes(fill = coverage), position = "stack", width = 1) +
    viridis::scale_fill_viridis(na.value = "grey") +
    facet_grid(celltype_label ~ .,
               scales = ifelse(y_var == "percent_cells", "fixed", "free_y")) +
    theme_minimal() +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.direction = "horizontal") +
    labs(title = p_title, subtitle = p_subtitle, y = y_var) +
    guides(colour = "none")
  if (facet_by_region) {
    p <-
      p +
      facet_grid(celltype_label ~ region, space = "free_x",
                  scales = ifelse(y_var == "percent_cells", "free_x", "free")) +
      scale_x_discrete(expand = c(0, 0)) +
      p_facet_theme
  }
  p
}

plot_ge <- function(p_dat, facet_by_region = FALSE) {
  p <-
    p_dat %>%
    dplyr::summarise(`% of cells\nw/ cov` = 100 *
                        sum(n_cells[coverage >= 2]) / sum(n_cells)) %>%
    ggplot(aes(x = pos, y = 1, fill = `% of cells\nw/ cov`,
                colour = `% of cells\nw/ cov`)) +
    geom_tile(width = 1, linewidth = 0.001) +
    viridis::scale_fill_viridis() +
    viridis::scale_colour_viridis() +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.direction = "horizontal")
  if (facet_by_region) {
    p <-
      p +
      facet_grid(. ~ region, scales = "free", space = "free_x") +
      scale_x_discrete(expand = c(0, 0)) +
      p_facet_theme
  }
  p
}

plot_mut_pie <- function(cov_per_ct, mut_annots, g, p_title) {

  # sort ct order
  ct_order <-
    mut_annots %>%
    dplyr::filter(!is.na(celltype)) %>%
    dplyr::pull(celltype) %>%
    {c("unannotated", ., "all")} %>%
    unique()

  # get mut sites only
  mut_sites <-
    mut_annots %>%
    dplyr::filter(!is.na(mut_ref), run == opts$meta_run) %>%
    dplyr::mutate(
      region = mut_ref,
      radius = 0.45,
      celltype = factor(celltype, levels = ct_order),
      celltype_n = as.integer(celltype),
      arr = ifelse(g$strand == "+", pos, -pos),
      mut_ref = forcats::fct_reorder(mut_ref, arr),
      mut_ref_n = as.integer(mut_ref),
      label = paste0(ifelse(total_depth == 0, "-", round(100 * alt_depth / total_depth, 1)),
                     "%\n", alt_depth, "/", total_depth),
      `variant detected?` = ifelse(alt_depth > 0, "yes", NA)) %>%
    dplyr::select(-total_depth)

  # plot mut pie
  p_mut_pie <-
    mut_sites %>%
    ggplot(aes(x = mut_ref_n, y = celltype_n)) +
    geom_tile(aes(colour = `variant detected?`), fill = NA, linewidth = 1) +
    scatterpie::geom_scatterpie(
      aes(x = mut_ref_n, y = celltype_n, r = radius),
      cols = c("ref_depth", "alt_depth"),
      legend_name = "depth", na.rm = TRUE, colour = NA) +
    geom_text(aes(label = label)) +
    coord_equal() +
    scale_x_continuous(breaks = unique(mut_sites$mut_ref_n),
                       labels = unique(mut_sites$mut_ref),
                       guide = guide_axis(angle = -45),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = unique(mut_sites$celltype_n),
                       labels = unique(mut_sites$celltype),
                       expand = c(0, 0)) +
    theme_minimal() +
    scale_fill_manual(values = c("ref_depth" = "#c0d7fa",
                                 "alt_depth" = "#eabac7")) +
    scale_colour_manual(values = c("yes" = "red"), na.translate = FALSE,
                        na.value = NA) +
    labs(x = "variant", title = p_title,
          subtitle = "Proportion of alt reads / all reads at the site") +
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_line())

  # plot ge, by celltype
  p_mut_ge <-
    cov_per_ct %>%
    dplyr::mutate(celltype = factor(celltype, levels = ct_order)) %>%
    dplyr::inner_join(mut_sites, by = c("celltype", "chr", "pos", "gene")) %>%
    dplyr::group_by(mut_ref, celltype, `variant detected?`) %>%
    dplyr::summarise(
      n_cells_w_cov = sum(n_cells[coverage >= 1]),
      n_cells_total = sum(n_cells),
      `% of cells\nw/ reads` = 100 * n_cells_w_cov / n_cells_total) %>%
    ggplot(aes(x = mut_ref, y = celltype)) +
    geom_tile(aes(fill = `% of cells\nw/ reads`), width = 1) +
    geom_tile(aes(colour = `variant detected?`), fill = NA, linewidth = 1) +
    geom_label(aes(label = paste0(round(`% of cells\nw/ reads`, 1), "%\n",
                                  n_cells_w_cov, "/", n_cells_total))) +
    viridis::scale_fill_viridis() +
    scale_colour_manual(values = c("yes" = "red"), na.translate = FALSE,
                        na.value = NA) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_discrete(expand = c(0, 0), guide = guide_axis(angle = -45)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_equal() +
    labs(subtitle = "Proportion of cells with at least 1 read at the site")

  # plot distance from polyA
  p_dist_poly_a <-
    mut_sites %>%
    dplyr::mutate(pos_polyA = ifelse(g$strand == "+", g$end, g$start),
                  dist_polyA = abs(pos - pos_polyA),
                  lab_polyA = prettyNum(dist_polyA, big.mark = ",",
                                        scientific = FALSE)) %>%
    dplyr::distinct(mut_ref, dist_polyA, lab_polyA) %>%
    ggplot(aes(x = mut_ref, y = dist_polyA)) +
    geom_col() +
    geom_text(aes(label = lab_polyA), vjust = 1.5) +
    theme_void() +
    labs(subtitle = "Genomic distance from polyA site") +
    scale_y_reverse()

  # patchwork
  p_mut_pie / p_mut_ge / p_dist_poly_a
}

# initialise list of lists plots
p <- list()
p[["genic"]] <- p[["exonic"]] <- p[["ccds"]] <- list()
p_subtitle <- paste0(opts$meta_run, " ", opts$gene, " (", g$strand, " strand, ",
                     prettyNum(g$end - g$start, big.mark = ",",
                               scientific = FALSE),
                     "bp, min coverage = ", opts$min_cov, ")")

# plot mutations from Mitchell 2022
if (nrow(mut_annots) > 0) {
  p[["genic"]][["wgs_mutations"]] <-
    plot_wgs_mutations(
      mut_annots,
      muts_subtitle = "SNVs called in WGS from Mitchell 2022")
  p[["exonic"]][["wgs_mutations"]] <-
    plot_wgs_mutations(
      dplyr::right_join(mut_annots, regions[["exonic"]]),
      facet_by_region = TRUE,
      muts_subtitle = "exonic SNVs called in WGS from Mitchell 2022")
  p[["ccds"]][["wgs_mutations"]] <-
    plot_wgs_mutations(
      dplyr::right_join(mut_annots, regions[["ccds"]]),
      facet_by_region = TRUE,
      muts_subtitle = "cCDS SNVs called in WGS from Mitchell 2022")
  p_wgs_muts_height <- 0.5
} else {
  p[["genic"]][["wgs_mutations"]] <- ggplot() + theme_void()
  p[["exonic"]][["wgs_mutations"]] <- ggplot() + theme_void()
  p[["ccds"]][["wgs_mutations"]] <- ggplot() + theme_void()
  p_wgs_muts_height <- 0.01
}

# plot mutations from PB panel
if (opts$meta_kit == "PB" & opts$meta_seq_type == "panel" & nrow(mut_annots) > 0) {

  # plot mutations
  p[["genic"]][["mutations"]] <-
    plot_mutations(mut_annots,
                   muts_subtitle = "mutations called")
  p[["exonic"]][["mutations"]] <-
    plot_mutations(dplyr::right_join(mut_annots, regions[["exonic"]]),
                   facet_by_region = TRUE,
                   muts_subtitle = "exonic mutations called")
  p[["ccds"]][["mutations"]] <-
    plot_mutations(dplyr::right_join(mut_annots, regions[["ccds"]]),
                   facet_by_region = TRUE,
                   muts_subtitle = "cCDS mutations called")
  p_muts_height <- 0.5

  # plot mut pie, save
  p_mut_pie <- plot_mut_pie(cov_per_ct, mut_annots, g, p_subtitle)
  saveRDS(p_mut_pie, paste0(opts$meta_run, "_", opts$gene,
                            "_pie_mutations_plot.rds"))
  ggsave(
    paste0(opts$meta_run, "_", opts$gene, "_pie_mutations_plot.",
           opts$plot_device),
    p_mut_pie,
    dpi = 500, width = 18, height = 18)

} else {
  p[["genic"]][["mutations"]] <- ggplot() + theme_void()
  p[["exonic"]][["mutations"]] <- ggplot() + theme_void()
  p[["ccds"]][["mutations"]] <- ggplot() + theme_void()
  p_muts_height <- 0.01
}

# plot transcripts
p[["genic"]][["transcripts"]] <- plot_transcripts(features)
p[["exonic"]][["transcripts"]] <-
  features %>%
  dplyr::group_by(dplyr::across(-c(start, end))) %>%
  dplyr::reframe(pos = start:end) %>%
  dplyr::right_join(regions[["exonic"]] %>% dplyr::select(chr, region, pos)) %>%
  dplyr::group_by(dplyr::across(-c(pos))) %>%
  dplyr::summarise(start = min(pos), end = max(pos)) %>%
  plot_transcripts(facet_by_region = TRUE)
p[["ccds"]][["transcripts"]] <-
  features %>%
  dplyr::group_by(dplyr::across(-c(start, end))) %>%
  dplyr::reframe(pos = start:end) %>%
  dplyr::right_join(regions[["ccds"]]) %>%
  dplyr::group_by(dplyr::across(-c(pos))) %>%
  dplyr::summarise(start = min(pos), end = max(pos)) %>%
  plot_transcripts(facet_by_region = TRUE)

# plot coverage, n cells
p[["genic"]][["cov"]] <-
  plot_cov(cov_per_ct)
p[["exonic"]][["cov"]] <-
  plot_cov(dplyr::right_join(cov_per_ct, regions[["exonic"]]),
           facet_by_region = TRUE, p_title = "Exonic coverage per cell")
p[["ccds"]][["cov"]] <-
  plot_cov(dplyr::right_join(cov_per_ct, regions[["ccds"]]),
           facet_by_region = TRUE, p_title = "Canonical CDS coverage per cell")

# plot coverage, % cells
p[["genic"]][["cov_pct"]] <-
  plot_cov(cov_per_ct, y_var = "percent_cells")
p[["exonic"]][["cov_pct"]] <-
  plot_cov(dplyr::right_join(cov_per_ct, regions[["exonic"]]),
           facet_by_region = TRUE, y_var = "percent_cells",
           p_title = "Exonic coverage per cell")
p[["ccds"]][["cov_pct"]] <-
  plot_cov(dplyr::right_join(cov_per_ct, regions[["ccds"]]),
           facet_by_region = TRUE, y_var = "percent_cells",
           p_title = "Canonical CDS coverage per cell")

# plot genotyping efficiency
p[["genic"]][["ge"]] <-
  cov_per_ct %>%
  dplyr::group_by(pos) %>%
  plot_ge()
p[["exonic"]][["ge"]] <-
  cov_per_ct %>%
  dplyr::inner_join(regions[["exonic"]]) %>%
  dplyr::group_by(region, pos) %>%
  plot_ge(facet_by_region = TRUE)
p[["ccds"]][["ge"]] <-
  cov_per_ct %>%
  dplyr::inner_join(regions[["ccds"]]) %>%
  dplyr::group_by(region, pos) %>%
  plot_ge(facet_by_region = TRUE)

# if strand is "-", flip the x axis in all plots
if (g$strand == "-") {
  p <-
    purrr::map(p, function(p1) {
      purrr::map(p1, function(p2) {
        p2 + scale_x_reverse(expand = c(0, 0))
      })
    })
}

# save plots
purrr::walk2(names(p), p, function(lvl, p1) {
  ggsave(
    paste0(opts$meta_run, "_", opts$gene, "_", lvl, "_cov_plot.",
           opts$plot_device),
    p1[["wgs_mutations"]] /
    p1[["mutations"]] /
    p1[["cov"]] /
    p1[["ge"]] /
    p1[["transcripts"]] +
    plot_layout(heights = c(p_wgs_muts_height, p_muts_height, 4, 0.5, 2)),
    dpi = 500)
  ggsave(
    paste0(opts$meta_run, "_", opts$gene, "_", lvl, "_cov_pct_plot.",
           opts$plot_device),
    p1[["wgs_mutations"]] /
    p1[["mutations"]] /
    p1[["cov_pct"]] /
    p1[["ge"]] /
    p1[["transcripts"]] +
    plot_layout(heights = c(p_wgs_muts_height, p_muts_height, 4, 0.5, 2)),
    dpi = 500)
})