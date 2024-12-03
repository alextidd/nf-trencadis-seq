#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(ggplot2)
library(patchwork)
library(optparse)

# parse arguments
option_list <- list(
  make_option("--meta_id", type = "character"),
  make_option("--gene", type = "character"),
  make_option("--gene_bed", type = "character"),
  make_option("--gene_regions", type = "character"),
  make_option("--gene_features", type = "character"),
  make_option("--cov_per_ct", type = "character"),
  make_option("--geno_per_ct", type = "character"),
  make_option("--min_cov", type = "numeric"),
  make_option("--plot_device", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")
print(opts)

# read in genotyped mutations (must explicitly set col_types because ref/alt 
# values of T will be interpreted as boolean column of TRUE)
geno_per_ct <-
  readr::read_tsv(opts$geno_per_ct,
                  col_types = readr::cols(ref = readr::col_character(),
                                          alt = readr::col_character()))

# read coverage data
cov_per_ct <- readRDS(opts$cov_per_ct)

# read gene data
g <- readr::read_tsv(opts$gene_bed,
                     col_names = c("chr", "start", "end", "gene", "strand"))
gene_regions <- readRDS(opts$gene_regions)
gene_features <- readRDS(opts$gene_features)

# functions
p_facet_theme <-
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "grey", fill = NA,
                                    linewidth = 0.1),
        strip.background = element_rect(color = "grey", fill = NA,
                                        linewidth = 0.1, linetype = "solid"))

plot_mutations <-
  function(p_dat, facet_by_region = FALSE, facet_by_celltype = FALSE,
           muts_subtitle = "mutations") {

    mut_type_pal <- c("snv" = "#241623", "ins" = "#D0CD94", "del" = "#68B0B6")
    n_muts <- dplyr::n_distinct(p_dat$mut_id, na.rm = TRUE)
    p <-
      p_dat %>%
      dplyr::select(dplyr::any_of(c("pos", "mut_type", "ref", "alt",
                                    "region", "celltype"))) %>%
      dplyr::mutate(mut = paste(ref, alt, sep = ">")) %>%
      dplyr::distinct() %>%
      ggplot(aes(x = pos, y = 0, colour = mut_type)) +
      geom_hline(yintercept = 0, colour = "grey") +
      geom_point(na.rm = TRUE) +
      theme_void() +
      scale_colour_manual(values = mut_type_pal,
                          na.translate = FALSE, na.value = NA) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      labs(subtitle = paste(n_muts, muts_subtitle))

    if (facet_by_region & facet_by_celltype) {
      p <-
        p +
        facet_grid(celltype ~ region, scales = "free_x", space = "free_x") +
        p_facet_theme
    } else if (facet_by_region & !facet_by_celltype) {
      p <-
        p +
        facet_grid(. ~ region, scales = "free_x", space = "free_x") +
        p_facet_theme
    } else if (!facet_by_region & facet_by_celltype) {
      p <-
        p +
        facet_grid(celltype ~ .) +
        p_facet_theme
    }

    p
  }

plot_genotyped_mutations <-
  function(p_dat, facet_by_region = FALSE, facet_by_celltype = FALSE,
           muts_subtitle = "mutations") {
    p_dat %>%
      dplyr::mutate(
        mut_type = dplyr::case_when(alt_depth == 0 | is.na(alt_depth) ~ NA,
                                    TRUE ~ mut_type),
        mut_id = dplyr::case_when(alt_depth == 0 | is.na(alt_depth) ~ NA,
                                  TRUE ~ mut_id)) %>%
      plot_mutations(facet_by_region = facet_by_region,
                     facet_by_celltype = facet_by_celltype,
                     muts_subtitle = muts_subtitle)
  }

plot_transcripts <- function(p_dat, facet_by_region = FALSE) {
  p <-
    p_dat %>%
    dplyr::filter(type != "gene") %>%
    ggplot(aes(x = start, xend = end,
               y = transcript_id, yend = transcript_id, linewidth = type)) +
    geom_segment(data = . %>% dplyr::filter(type %in% c("transcript", "exon")),
                  aes(colour = exon_alternating)) +
    scale_linewidth_manual(values = c("transcript" = 0.5, "exon" = 3)) +
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

plot_mut_pie <- function(cov_per_ct, geno_per_ct, g, p_title) {

  # sort ct order
  ct_order <-
    geno_per_ct %>%
    dplyr::filter(!is.na(celltype)) %>%
    dplyr::pull(celltype) %>%
    {c("unannotated", ., "all")} %>%
    unique()

  # get mut sites only
  mut_sites <-
    geno_per_ct %>%
    dplyr::filter(!is.na(mut_id)) %>%
    dplyr::mutate(
      region = mut_id,
      radius = 0.45,
      celltype = factor(celltype, levels = ct_order),
      celltype_n = as.integer(celltype),
      arr = ifelse(g$strand == "+", pos, -pos),
      mut_id = forcats::fct_reorder(mut_id, arr),
      mut_ref_n = as.integer(mut_id),
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
                       labels = unique(mut_sites$mut_id),
                       guide = guide_axis(angle = -45),
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = unique(mut_sites$celltype_n),
                       labels = unique(mut_sites$celltype),
                       expand = c(0.01, 0.01)) +
    theme_minimal() +
    scale_fill_manual(values = c("ref_depth" = "#c0d7fa",
                                 "alt_depth" = "#eabac7")) +
    scale_colour_manual(values = c("yes" = "#b30000"), na.translate = FALSE,
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
    dplyr::group_by(mut_id, celltype, `variant detected?`) %>%
    dplyr::summarise(
      n_cells_w_cov = sum(n_cells[coverage >= 1]),
      n_cells_total = sum(n_cells),
      `% of cells\nw/ reads` = 100 * n_cells_w_cov / n_cells_total) %>%
    ggplot(aes(x = mut_id, y = celltype)) +
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

  # patchwork
  p_mut_pie / p_mut_ge
}

# initialise list of lists plots
p <- list()
p[["genic"]] <- p[["exonic"]] <- p[["ccds"]] <- list()
p_subtitle <- paste0(opts$meta_id, " ", opts$gene, " (", g$strand, " strand, ",
                     prettyNum(g$end - g$start, big.mark = ",",
                               scientific = FALSE),
                     "bp, min coverage = ", opts$min_cov, ")")

# plot mutations from Mitchell 2022
if (nrow(geno_per_ct) > 0) {
  p[["genic"]][["wgs_mutations"]] <-
    plot_mutations(
      dplyr::right_join(geno_per_ct, gene_regions[["genic"]]),
      muts_subtitle = "mutations")
  p[["exonic"]][["wgs_mutations"]] <-
    plot_mutations(
      dplyr::right_join(geno_per_ct, gene_regions[["exonic"]]),
      facet_by_region = TRUE,
      muts_subtitle = "exonic mutations")
  p[["ccds"]][["wgs_mutations"]] <-
    plot_mutations(
      dplyr::right_join(geno_per_ct, gene_regions[["ccds"]]),
      facet_by_region = TRUE,
      muts_subtitle = "cCDS mutations")
  p_wgs_muts_height <- 0.5
} else {
  p[["genic"]][["wgs_mutations"]] <- ggplot() + theme_void()
  p[["exonic"]][["wgs_mutations"]] <- ggplot() + theme_void()
  p[["ccds"]][["wgs_mutations"]] <- ggplot() + theme_void()
  p_wgs_muts_height <- 0.01
}

# plot mutations from PB panel
if (sum(geno_per_ct$alt_depth) > 0) {

  # plot mutations
  p[["genic"]][["mutations"]] <-
    plot_genotyped_mutations(
      dplyr::right_join(geno_per_ct, gene_regions[["genic"]]),
      facet_by_celltype = TRUE,
      muts_subtitle = "mutations genotyped")
  p[["exonic"]][["mutations"]] <-
    plot_genotyped_mutations(
      dplyr::right_join(geno_per_ct, gene_regions[["exonic"]]),
      facet_by_region = TRUE, facet_by_celltype = TRUE,
      muts_subtitle = "exonic mutations genotyped")
  p[["ccds"]][["mutations"]] <-
    plot_genotyped_mutations(
      dplyr::right_join(geno_per_ct, gene_regions[["ccds"]]),
      facet_by_region = TRUE, facet_by_celltype = TRUE,
      muts_subtitle = "cCDS mutations genotyped")
  p_muts_height <- 1.5

  # plot mut pie, save
  p_mut_pie <- plot_mut_pie(cov_per_ct, geno_per_ct, g, p_subtitle)
  saveRDS(p_mut_pie, paste0(opts$meta_id, "_", opts$gene,
                            "_pie_mutations_plot.rds"))
  ggsave(
    paste0(opts$meta_id, "_", opts$gene, "_pie_mutations_plot.",
           opts$plot_device),
    p_mut_pie, dpi = 500)

} else {
  p[["genic"]][["mutations"]] <- ggplot() + theme_void()
  p[["exonic"]][["mutations"]] <- ggplot() + theme_void()
  p[["ccds"]][["mutations"]] <- ggplot() + theme_void()
  p_muts_height <- 0.01
}

# plot transcripts
p[["genic"]][["transcripts"]] <- plot_transcripts(gene_features)
p[["exonic"]][["transcripts"]] <-
  gene_features %>%
  dplyr::group_by(dplyr::across(-c(start, end))) %>%
  dplyr::reframe(pos = start:end) %>%
  dplyr::right_join(gene_regions[["exonic"]] %>%
                      dplyr::select(chr, region, pos)) %>%
  dplyr::group_by(dplyr::across(-c(pos))) %>%
  dplyr::summarise(start = min(pos), end = max(pos)) %>%
  plot_transcripts(facet_by_region = TRUE)
p[["ccds"]][["transcripts"]] <-
  gene_features %>%
  dplyr::group_by(dplyr::across(-c(start, end))) %>%
  dplyr::reframe(pos = start:end) %>%
  dplyr::right_join(gene_regions[["ccds"]]) %>%
  dplyr::group_by(dplyr::across(-c(pos))) %>%
  dplyr::summarise(start = min(pos), end = max(pos)) %>%
  plot_transcripts(facet_by_region = TRUE)

# plot coverage, n cells
p[["genic"]][["cov"]] <-
  plot_cov(cov_per_ct)
p[["exonic"]][["cov"]] <-
  plot_cov(dplyr::right_join(cov_per_ct, gene_regions[["exonic"]]),
           facet_by_region = TRUE, p_title = "Exonic coverage per cell")
p[["ccds"]][["cov"]] <-
  plot_cov(dplyr::right_join(cov_per_ct, gene_regions[["ccds"]]),
           facet_by_region = TRUE, p_title = "Canonical CDS coverage per cell")

# plot coverage, % cells
p[["genic"]][["cov_pct"]] <-
  plot_cov(cov_per_ct, y_var = "percent_cells")
p[["exonic"]][["cov_pct"]] <-
  plot_cov(dplyr::right_join(cov_per_ct, gene_regions[["exonic"]]),
           facet_by_region = TRUE, y_var = "percent_cells",
           p_title = "Exonic coverage per cell")
p[["ccds"]][["cov_pct"]] <-
  plot_cov(dplyr::right_join(cov_per_ct, gene_regions[["ccds"]]),
           facet_by_region = TRUE, y_var = "percent_cells",
           p_title = "Canonical CDS coverage per cell")

# plot genotyping efficiency
p[["genic"]][["ge"]] <-
  cov_per_ct %>%
  dplyr::group_by(pos) %>%
  plot_ge()
p[["exonic"]][["ge"]] <-
  cov_per_ct %>%
  dplyr::inner_join(gene_regions[["exonic"]]) %>%
  dplyr::group_by(region, pos) %>%
  plot_ge(facet_by_region = TRUE)
p[["ccds"]][["ge"]] <-
  cov_per_ct %>%
  dplyr::inner_join(gene_regions[["ccds"]]) %>%
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
    paste0(opts$meta_id, "_", opts$gene, "_", lvl, "_cov_plot.",
           opts$plot_device),
    p1[["wgs_mutations"]] /
    p1[["mutations"]] /
    p1[["cov"]] /
    p1[["ge"]] /
    p1[["transcripts"]] +
    plot_layout(heights = c(p_wgs_muts_height, p_muts_height, 4, 0.5, 2)),
    dpi = 500)
  ggsave(
    paste0(opts$meta_id, "_", opts$gene, "_", lvl, "_cov_pct_plot.",
           opts$plot_device),
    p1[["wgs_mutations"]] /
    p1[["mutations"]] /
    p1[["cov_pct"]] /
    p1[["ge"]] /
    p1[["transcripts"]] +
    plot_layout(heights = c(p_wgs_muts_height, p_muts_height, 4, 0.5, 2)),
    dpi = 500)
})
