#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# parse arguments
option_list <- list(
  make_option("--gene", type = "character"),
  make_option("--gene_bed", type = "character"),
  make_option("--exonic_bed", type = "character"),
  make_option("--refcds", type = "character"),
  make_option("--gene_features", type = "character"),
  make_option("--cov", type = "character"),
  make_option("--cell_barcodes", type = "character"),
  make_option("--celltypes", type = "character"),
  make_option("--meta_run", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# print the options
print(opts)

# get gene info (incl. strand to orient 5'/3' ends)
g <- readr::read_tsv(opts$gene_bed,
                     col_names = c("chr", "start", "end", "gene", "strand"))

# create list of regions
regions <- list()

# get exonic regions
regions[["exonic"]] <-
  readr::read_tsv(opts$exonic_bed,
                  col_names = c("chr", "start", "end")) %>%
  # label each region, expand to all positions
  dplyr::distinct() %>%
  dplyr::mutate(arr = dplyr::case_when(g$strand == "+" ~ start,
                                        g$strand == "-" ~ -start)) %>%
  dplyr::arrange(arr) %>%
  dplyr::mutate(region = dplyr::row_number()) %>%
  dplyr::group_by(chr, region) %>%
  dplyr::reframe(pos = start:end) %>%
  # get exonic distance from 3'
  dplyr::group_by(pos) %>%
  dplyr::mutate(exonic_distance_from_5_prime = dplyr::cur_group_id()) %>%
  dplyr::ungroup()

# get ccds regions
load(opts$refcds)
ref_gene <-
  RefCDS[[which(purrr::map_lgl(RefCDS, ~ .x$gene_name == opts$gene))]]
regions[["ccds"]] <-
  tibble::tibble(chr = paste0("chr", ref_gene$chr),
                 start = ref_gene$intervals_cds[, 1],
                 end = ref_gene$intervals_cds[, 2])  %>%
  dplyr::mutate(arr = dplyr::case_when(g$strand == "+" ~ start,
                                        g$strand == "-" ~ -start)) %>%
  dplyr::arrange(arr) %>%
  dplyr::mutate(region = dplyr::row_number()) %>%
  dplyr::group_by(chr, region) %>%
  dplyr::reframe(pos = start:end)

# get features (must wrangle from GFF3 attributes)
features <-
  readr::read_tsv(opts$gene_features) %>%
  dplyr::mutate(
    id = gsub(";.*", "", attributes)) %>%
  tidyr::separate_longer_delim(cols = "attributes", delim = ";") %>%
  tidyr::separate_wider_delim(
    cols = "attributes", delim = "=", names = c("name", "value")) %>%
  dplyr::filter(name %in% c("exon_id", "transcript_id", "transcript_type",
                            "exon_number")) %>%
  tidyr::pivot_wider() %>%
  readr::type_convert() %>%
  dplyr::mutate(
    # move protein-coding to top of plot
    transcript_type = forcats::fct_relevel(transcript_type, "protein_coding"),
    # alternate exons for colouring
    exon_alternating = as.character(exon_number %% 2))

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

# save
g %>% saveRDS(paste0(opts$meta_run, "_", opts$gene, ".rds"))
regions %>% saveRDS(paste0(opts$meta_run, "_", opts$gene, "_regions.rds"))
features %>% saveRDS(paste0(opts$meta_run, "_", opts$gene, "_features.rds"))
cov_per_ct %>% saveRDS(paste0(opts$meta_run, "_", opts$gene, "_cov_per_ct.rds"))
exonic_cov %>% readr::write_tsv(paste0(opts$meta_run, "_", opts$gene,
                                       "_coverage_per_exonic_position.tsv"))