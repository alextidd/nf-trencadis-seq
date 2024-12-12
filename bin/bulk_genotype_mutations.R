#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# parse arguments
option_list <- list(
  make_option("--bam", type = "character"),
  make_option("--gene_bed", type = "character"),
  make_option("--mutations", type = "character"),
  make_option("--min_MQ", type = "numeric"),
  make_option("--min_BQ", type = "numeric"))
opts <- parse_args(OptionParser(option_list = option_list))
saveRDS(opts, "opts.rds")
# opts <- readRDS("opts.rds")

# load the mutations in the gene
gene_pos <-
  readr::read_tsv(opts$gene_bed,
                  col_names = c("chr", "start", "end", "gene", "strand")) %>%
  dplyr::mutate(chr = as.character(chr)) %>%
  dplyr::group_by(chr, gene) %>%
  dplyr::reframe(pos = start:end)
mutations <-
  readr::read_tsv(opts$mutations) %>%
  dplyr::distinct(chr, pos, ref, alt) %>%
  dplyr::inner_join(gene_pos) %>%
  dplyr::mutate(
    mut_id = paste(chr, pos, ref, alt, sep = "-"),
    mut_type = dplyr::case_when(nchar(ref) == 1 & nchar(alt) == 1 ~ "snv",
                                nchar(ref) == 1 & nchar(alt) > 1 ~ "ins",
                                nchar(ref) > 1 & nchar(alt) == 1 ~ "del",
                                TRUE ~ "complex"))

# check for mutations that are not snv / ins / del                            
if ("complex" %in% mutations$mut_type) {
  message("Complex mutations are not supported!")
  message(paste(mutations %>% dplyr::filter(mut_type == "complex") %>% nrow(),
                "complex mutation(s) found and will be removed."))
  mutations <- mutations %>% dplyr::filter(mut_type != "complex")
} 

if (nrow(mutations) > 0) {

  geno <-
    mutations %>%
    dplyr::distinct(chr, pos, gene, ref, alt, mut_id, mut_type) %>%
    purrr::pmap(function(chr, pos, gene, ref, alt, mut_id, mut_type) {

      # print mut
      paste(chr, pos, ref, alt, mut_type, "\n") %>% cat()
      
      # query bam
      calls <- deepSNV::bam2R(opts$bam, chr, pos, pos,
                              q = opts$min_BQ, mq = opts$min_MQ)
      
      # count all reads at site
      total_depth <- sum(calls[, c("A", "C", "G", "T", "a", "c", "g", "t",
                                    "DEL", "INS", "del", "ins")],
                          na.rm = TRUE)
    
      # count ref reads at site
      # take first character (in case it is a deletion)
      ref_1 <- substr(ref, 1, 1)
      ref_depth <- sum(calls[, c(ref_1, tolower(ref_1))], na.rm = TRUE)

      # count mutant reads at site
      if (mut_type == "snv") {
        # count mutant reads at site
        alt_depth <- sum(calls[, c(alt, tolower(alt))], na.rm = TRUE)
      } else if (mut_type %in% c("ins", "del")) {
        # count ins or del reads at site (don't check sequence)
        alt_depth <- sum(calls[, c(mut_type, toupper(mut_type))], na.rm = TRUE)
      }
      
      tibble::tibble(chr = chr, pos = pos, ref = ref, alt = alt,
                      mut_id = mut_id, mut_type = mut_type, gene = gene,
                      total_depth = total_depth, ref_depth = ref_depth,
                      alt_depth = alt_depth,
                      alt_vaf = alt_depth / total_depth)
      
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(alt_vaf = alt_depth / total_depth)

} else {

  geno <- tibble::tibble()

}

# write
mutations %>% readr::write_tsv("mutations.tsv")
geno %>% readr::write_tsv("genotyped_mutations.tsv")