#!/usr/bin/env Rscript

# libraries
library(GenomicRanges)
library(magrittr)
library(optparse)

# parse arguments
option_list <- list(
  make_option("--gene", type = "character"),
  make_option("--gencode_gff3", type = "character"),
  make_option("--refcds", type = "character"))
opts <- parse_args(OptionParser(option_list = option_list))
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

print(opts)

print("get ccds coords from refcds")
load(opts$refcds)
g <- RefCDS[[which(purrr::map_lgl(RefCDS, ~ .x$gene_name == opts$gene))]]
ccds <-
  tibble::tibble(chr = paste0("chr", g$chr), start = g$intervals_cds[, 1],
                 end = g$intervals_cds[, 2], type = "ccds")

print("get feature coords from the gff")
gff <- ape::read.gff(opts$gencode_gff3)
features <-
  gff %>%
  dplyr::filter(
    grepl(paste0("gene_name=", opts$gene, ";"), attributes),
    # get protein-coding gene coords
    (type == "gene" & grepl("gene_type=protein_coding;", attributes)) |
    # get transcript / exon coords
    (type %in% c("transcript", "exon"))) %>%
  dplyr::transmute(chr = seqid, start, end, gene = opts$gene, strand, type,
                   attributes) %>%
  dplyr::bind_rows(ccds) %>%
  # extract info from attributes
  dplyr::mutate(
    id = gsub(";.*", "", attributes)) %>%
  tidyr::separate_longer_delim(cols = "attributes", delim = ";") %>%
  tidyr::separate_wider_delim(
    cols = "attributes", delim = "=", names = c("name", "value")) %>%
  dplyr::filter(name %in% c("gene_id", "exon_id", "transcript_id",
                            "transcript_type", "exon_number")) %>%
  tidyr::pivot_wider() %>%
  readr::type_convert() %>%
  dplyr::mutate(
    # move protein-coding to top of plot
    transcript_type = forcats::fct_relevel(transcript_type, "protein_coding"),
    # alternate exons for colouring
    exon_alternating = as.character(exon_number %% 2))

print("extract the gene coords from the gff")
gene_bed <-
  features %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(chr, start, end, gene, strand)  

# create list of regions
regions <- list()

print("get genic regions")
regions[["genic"]] <-
  gene_bed %>%
  dplyr::group_by(chr, gene, strand) %>%
  dplyr::reframe(pos = start:end)

print("get exonic regions")
regions[["exonic"]] <-
  features %>%
  # collapse all exonic regions
  dplyr::filter(type == "exon") %>%
  {GRanges(seqnames = .$chr,
            ranges = IRanges(.$start, .$end))} %>%
  reduce() %>%
  {tibble::tibble(chr = as.vector(seqnames(.)),
                  start = start(.),
                  end = end(.))} %>%
  # label each region, expand to all positions
  dplyr::distinct() %>%
  dplyr::mutate(arr = dplyr::case_when(gene_bed$strand == "+" ~ start,
                                       gene_bed$strand == "-" ~ -start)) %>%
  dplyr::arrange(arr) %>%
  dplyr::mutate(region = dplyr::row_number()) %>%
  dplyr::group_by(chr, region) %>%
  dplyr::reframe(pos = start:end) %>%
  # get exonic distance from 5'
  dplyr::group_by(pos) %>%
  dplyr::mutate(exonic_distance_from_5_prime = dplyr::cur_group_id()) %>%
  dplyr::ungroup()

print("get ccds regions")
regions[["ccds"]] <-
  tibble::tibble(chr = paste0("chr", g$chr),
                 start = g$intervals_cds[, 1],
                 end = g$intervals_cds[, 2])  %>%
  dplyr::mutate(arr = dplyr::case_when(gene_bed$strand == "+" ~ start,
                                       gene_bed$strand == "-" ~ -start)) %>%
  dplyr::arrange(arr) %>%
  dplyr::mutate(region = dplyr::row_number()) %>%
  dplyr::group_by(chr, region) %>%
  dplyr::reframe(pos = start:end)

print("save output")
readr::write_tsv(gene_bed, paste0(opts$gene, ".bed"), col_names = FALSE)
saveRDS(regions, paste0(opts$gene, "_regions.rds"))
saveRDS(features, paste0(opts$gene, "_features.rds"))
