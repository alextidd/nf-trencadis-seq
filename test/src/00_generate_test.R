# libraries
library(magrittr)
library(dplyr)

# dir
wd <- getwd()

# define test set
curr_id <- "PB_panel_KX004"

# get genotyping per cell
genos <-
  c("ASXL1", "DNMT3A", "TP53") %>%
  purrr::map(~
    file.path("/lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/nf-trencadis-seq/runs/",
              curr_id, "/", .x, "/genotyped_mutations_per_cell.tsv") %>%
    readr::read_tsv(col_types = "ccdcccccddddc")) %>%
  bind_rows()

# get 10 mutated, 10 unmutated, 10 no coverage, at least two genes
test_genos <-
  genos %>%
  mutate(detected = ifelse(alt_vaf > 0, "mutated", "unmutated")) %>%
  group_by(detected, gene) %>%
  slice_sample(n = 50) %>%
  ungroup()

# write to files
unique(test_genos$gene) %>% writeLines("test/data/genes.txt")
unique(test_genos$barcode) %>% writeLines("test/data/cell_barcodes.txt")
test_genos %>% distinct(chr, pos, ref, alt) %>% readr::write_tsv("test/data/mutations.tsv")
test_genos %>% select(barcode, celltype) %>% readr::write_csv("test/data/celltypes.csv")
tibble::tibble(id = curr_id,
               bam = file.path(wd, "test/data/PB_panel_KX004_subset.bam"),
               cell_barcodes = file.path(wd, "test/data/cell_barcodes.txt"),
               mutations = file.path(wd, "test/data/mutations.tsv"),
               celltypes = file.path(wd, "test/data/celltypes.csv")) %>%
  readr::write_csv("test/out/samplesheet.csv")
