#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# parse arguments
option_list <- list(
  make_option("--gene", type = "character"),
  make_option("--regions", type = "character"),
  make_option("--cov", type = "character"),
  make_option("--cell_barcodes", type = "character"),
  make_option("--celltypes", type = "character"),
  make_option("--meta_id", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# print the options
print(opts)

# get regions
regions <- readRDS(opts$regions)

# get coverage data
cov <- readr::read_tsv(opts$cov)

# get all cell barcodes, remove suffix
cell_barcodes <- readLines(opts$cell_barcodes) %>% gsub("-1$", "", .)

# get celltype annotations
celltypes <- readr::read_csv(opts$celltypes)

# get lists of celltype-annotated barcodes and list of all barcodes
barcodes_ls <- split(celltypes$barcode, celltypes$celltype)
barcodes_ls[["all"]] <- cell_barcodes

# get the counts of each unique value in each matrix row, per celltype
cov_per_ct <-
  barcodes_ls %>%
  purrr::map(function(CBs) {

    # get a matrix of reads in cells
    mat <-
      cov %>%
      dplyr::select(dplyr::any_of(CBs)) %>%
      as.matrix()

    # add columns for missing CBs (cells with 0 coverage)
    CBs0 <- setdiff(CBs, colnames(mat))
    mat0 <- matrix(0, nrow = nrow(mat), ncol = length(CBs0))
    colnames(mat0) <- CBs0
    mat <- cbind(mat, mat0)

    # get counts of each unique value in each row
    counts <-
      lapply(1:nrow(mat), function(x) {
        table(mat[x, ]) %>%
          tibble::enframe(name = "coverage", value = "n_cells") %>%
          cbind(cov[x, c("chr", "pos", "gene")], .)
      }) %>%
      dplyr::bind_rows() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(coverage = as.numeric(coverage),
                    n_cells = as.numeric(n_cells),
                    total_cells = length(CBs))

    # return
    counts

  }) %>%
  dplyr::bind_rows(.id = "celltype") %>%
  dplyr::mutate(celltype_label = paste0(celltype, " (", total_cells, ")") %>%
                  forcats::fct_reorder(total_cells, .desc = TRUE))

# save coverage
cov_per_ct %>% readr::write_tsv(paste0(opts$meta_id, "_", opts$gene,
                                       "_coverage_per_celltype.tsv"))

# get total n cells (across celltypes)
total_cells_all_cts <-
  cov_per_ct %>%
  dplyr::distinct(celltype_label, total_cells) %>%
  {sum(.$total_cells)}

# pile up coverage across celltypes at each exonic position
exonic_cov <-
  cov_per_ct %>%
  dplyr::inner_join(regions[["exonic"]]) %>%
  # get number of cells per coverage per position
  dplyr::group_by(chr, pos, gene, region, exonic_distance_from_5_prime,
                  coverage) %>%
  dplyr::summarise(n_cells = sum(n_cells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cells = total_cells_all_cts)

# pile up coverage across celltypes at each ccds position
ccds_cov <-
  cov_per_ct %>%
  dplyr::inner_join(regions[["ccds"]]) %>%
  # get number of cells per coverage per position
  dplyr::group_by(chr, pos, gene, region, coverage) %>%
  dplyr::summarise(n_cells = sum(n_cells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cells = total_cells_all_cts)

# save subsetted coverage
exonic_cov %>% readr::write_tsv(paste0(opts$meta_id, "_", opts$gene,
                                       "_coverage_per_exonic_position.tsv"))
ccds_cov %>% readr::write_tsv(paste0(opts$meta_id, "_", opts$gene,
                                     "_coverage_per_ccds_position.tsv"))