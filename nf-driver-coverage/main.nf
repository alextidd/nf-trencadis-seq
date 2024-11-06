#!/usr/bin/env nextflow
// TODO: fix called mutation annotations to fit + facet by celltype
// TODO: concat PB WT and panels together to get all reads together and recall mutations
// TODO: get Rashesh to try venn diagram of barcode propagation as a sanity check
// TODO: extract barcodes from the BAMs directly + recreate revcomp conversion to check where it fails
// TODO: add gene expression distribution alongside mutations / coverage plots
// TODO: create summary plot per gene
// TODO: do PB panel / WT UMIs versus TX WT UMIs per gene - highlight drivers in red
// TODO: add celltype marker genes to the expression dotplot (ask Laura)
// TODO: run collectHSMetrics on PB panel data to characterise bait capture efficiency
// TODO: check over the transcript IDs that Carol from pipeline sent - begin deconvolution in the dataset
// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// validate input parameters
if (params.validate_params) {
  validateParameters()
}

// Download a given sample's BAM from iRODS
// Then either retrieve the BAI or make one via indexing
// The maxForks of 10 was set after asking jc18 about best iRODS practices
process irods {
  tag "${meta.run}"
  maxForks 10
  label 'normal4core'

  input:
  tuple val(meta), val(bam),
        path(celltypes), path(mut_annots), path(cell_barcodes)
  
  output:
  tuple val(meta), path("${meta.run}.bam"), path("${meta.run}.bam.bai"),
        path(celltypes), path(mut_annots), path(cell_barcodes)
  
  script:
  """
  iget -K ${bam} ${meta.run}.bam
  if [[ `ils ${bam}.bai | wc -l` == 1 ]]
  then
      iget -K ${bam}.bai ${meta.run}.bam.bai
  else
      samtools index -@ ${task.cpus} ${meta.run}.bam
  fi
  """
}

// The equivalent of an irods download, but for a local copy of samplesheet
// Symlink the BAM/BAI appropriately so they're named the right thing for downstream
process local {
  tag "${meta.run}"
  maxForks 10
  label 'normal4core'
  errorStrategy = {task.attempt <= 1 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), val(bam),
        path(celltypes), path(mut_annots), path(cell_barcodes)

  output:
  tuple val(meta), path("${meta.run}.bam"), path("${meta.run}.bam.bai"),
        path(celltypes), path(mut_annots), path(cell_barcodes)

  script:
  """
  # create local symbolic link 
  ln -s ${bam} ${meta.run}.bam
  if [ -f "${bam}.bai" ] ; then
      ln -s ${bam}.bai ${meta.run}.bam.bai
  else
      samtools index -@ ${task.cpus} ${meta.run}.bam
  fi
  """
}

// check bam is not truncated before proceeding
process check_bam {
  tag "${meta.run}"
  label 'normal'
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mut_annots),
        path(cell_barcodes)

  output:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mut_annots),
        path(cell_barcodes)

  script:
  """
  module load samtools-1.19/python-3.12.0 
  samtools quickcheck $bam
  """
}

// get driver gene canonical cds
process get_driver_gene_ccds {
  tag "${gene}"
  label 'normal'
  publishDir "${params.out_dir}/genes/${gene}/",
    mode: "copy",
    pattern: "*_ccds_regions.bed"

  input:
  val(gene)
  path(refcds)

  output:
  tuple val(gene), path("${gene}_ccds_regions.bed")

  script:
  """
  #!/usr/bin/env Rscript

  # libraries
  library(magrittr)

  # get canonical cds for each gene
  load("${refcds}")
  g <- RefCDS[[which(purrr::map_lgl(RefCDS, ~ .x\$gene_name == "${gene}"))]]
  tibble::tibble(chr = paste0("chr", g\$chr), start = g\$intervals_cds[, 1],
                 end = g\$intervals_cds[, 2]) %>%
    readr::write_tsv("${gene}_ccds_regions.bed")
  """
}

// get driver gene coords
process get_driver_gene_coords {
  tag "${gene}"
  label 'normal'
  publishDir "${params.out_dir}/genes/${gene}/",
    mode: "copy",
    pattern: "*.bed"

  input:
  tuple val(gene), path(ccds_bed)
  path(gencode_gff3)

  output:
  tuple val(gene),
        path(ccds_bed),
        path("${gene}.bed"),
        path("${gene}_features.bed"),
        path("${gene}_exonic_regions.bed"), emit: gene_features
  path "${gene}.bed", emit: gene_bed

  script:
  """
  # get the protein-coding gene coords from the gff3 file (chr start end gene strand)
  zcat ${gencode_gff3} \\
  | grep -w ${gene} \\
  | grep -w gene_type=protein_coding \\
  | awk -F"\\t" -v OFS="\\t" '\$3 == "gene" {print \$1, \$4, \$5, "${gene}", \$7}' \\
  > ${gene}.bed

  # get the transcript/exon coords from the gff3 file
  echo -e "chr\\tstart\\tend\\tgene\\ttype\\tattributes" > ${gene}_features.bed
  zcat ${gencode_gff3} \\
  | grep -w ${gene} \\
  | awk -F"\t" -v OFS="\t" '\$3 == "exon" || \$3 == "transcript" {print \$1, \$4, \$5, "${gene}", \$3, \$9}' \\
  >> ${gene}_features.bed

  # get all exonic regions, collapsing overlaps
  module load bedtools2-2.29.0/python-3.10.10
  cat ${gene}_features.bed \
  | awk -F"\\t" -v OFS="\\t" '\$5 == "exon" {print \$1, \$2, \$3}' \
  | sort -k1,1 -k2,2n \\
  | bedtools merge -i - \\
  > ${gene}_exonic_regions.bed
  """
}

// subset the BAMs to only the driver regions + only barcodes in the seurat dir
process subset_bams_to_drivers_and_barcodes {
  tag "${meta.run}_${gene}"
  label 'normal'
  publishDir "${params.out_dir}/runs/${meta.run}/${gene}/",
    mode: "copy",
    pattern: "*_subset.{bam,bam.bai}"

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mut_annots),
        path(cell_barcodes),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)

  output:
  tuple val(meta),
        path("${meta.run}_subset.bam"), path("${meta.run}_subset.bam.bai"),
        path(celltypes),
        path(mut_annots),
        path(cell_barcodes),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)

  script:
  """
  module load samtools-1.19/python-3.12.0 

  # subset to reads that fall in driver regions
  samtools view -L $gene_bed -b $bam > ${meta.run}_subset_drivers.bam
  samtools index -@ ${task.cpus} ${meta.run}_subset_drivers.bam

  # subset to cell barcodes
  subset-bam \
    -b ${meta.run}_subset_drivers.bam \
    --cell-barcodes $cell_barcodes \
    --out-bam ${meta.run}_subset.bam
  
  # index
  samtools index -@ ${task.cpus} ${meta.run}_subset.bam
  """
}

// split subsetted BAMs by cell, count coverage
process get_coverage_per_cell {
  tag "${meta.run}_${gene}"
  label 'week10gb'
  publishDir "${params.out_dir}/runs/${meta.run}/${gene}/",
    mode: "copy",
    pattern: "*_coverage_per_cell.tsv"

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mut_annots),
        path(cell_barcodes),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)

  output:
  tuple val(meta), path("${meta.run}_${gene}_coverage_per_cell.tsv"),
        path(celltypes),
        path(mut_annots),
        path(cell_barcodes),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)

  script:
  """
  module load samtools-1.19/python-3.12.0

  echo "getting cell barcodes"
  # get unique cell barcodes in the BAM, trim CB:Z: prefix
  samtools view $bam | cut -f12- | tr "\\t" "\\n" |
  grep "CB:Z:" | sed 's/^CB:Z://g' | awk '!x[\$0]++' \
  > cell_barcodes.txt

  echo "creating a bam per cell"
  rm -rf cell_bams ; mkdir -p cell_bams
  while read -r CB ; do
    # subset
    echo \$CB > \$CB.txt
    subset-bam \
      -b $bam \
      --cell-barcodes \$CB.txt \
      --out-bam cell_bams/\$CB.bam
    # index
    samtools index -@ ${task.cpus} cell_bams/\$CB.bam
    rm \$CB.txt
  done < cell_barcodes.txt

  # run per gene
  while read -r chr start end gene strand ; do

    echo "processing gene \$gene"

    # initialise gene coverage file
    rm -rf \${gene} ; mkdir \${gene}
    touch \${gene}_cov.tsv

    echo "calculating coverage per base per cell"
    while read -r CB ; do
      (
        echo "\$CB" ;
        samtools depth -a -r \$chr:\$start-\$end cell_bams/\$CB.bam \
        | cut -f3 ;
      ) > \${gene}/\$CB.tsv
    done < cell_barcodes.txt

    echo "creating coverage file"
    (
      echo -e "chr\\tpos\\tgene" ;
      echo -e "\$chr\\t\$start\\t\$end" \
      | awk -F"\\t" -v gene=\$gene -v OFS="\\t" \
        '{for(i=\$2; i<=\$3; i++) print \$1, i, gene}' ;
    ) > \${gene}_coords.tsv

    echo "adding coverage per cell to coverage file"
    # (sort for consistent ordering)
    paste \$(ls \${gene}/*.tsv | sort) > \${gene}_cov.tsv

    # combine coords and coverage
    paste \${gene}_coords.tsv \${gene}_cov.tsv > \${gene}_cov_per_cell.tsv

  done < $gene_bed

  # clean up cell bams dir
  rm -rf cell_bams/

  echo "combining outputs"
  ls *_cov_per_cell.tsv | head -1 | xargs head -1 \
  > ${meta.run}_${gene}_coverage_per_cell.tsv
  cat *_cov_per_cell.tsv | grep -vP "^chr\\tpos" \
  >> ${meta.run}_${gene}_coverage_per_cell.tsv
  """
}

// wrangle data for plots
process wrangle_data {
  tag "${meta.run}_${gene}"
  label 'normal20gb'
  maxRetries 10
  beforeScript "module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=~/R-tmp-4.4"
  publishDir "${params.out_dir}/genes/${gene}/",
    mode: "copy",
    pattern: "*_avg_coverage_per_ccds_base_per_cell.tsv"
  
  input:
  tuple val(meta), path(cov),
        path(celltypes),
        path(mut_annots),
        path(cell_barcodes),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)
  
  output:
  tuple val(meta), val(gene), path("*.rds"),
        emit: plot_data
  tuple val(gene),
        path("${meta.run}_${gene}_coverage_per_exonic_position.tsv"),
        emit: exonic_coverage
  path("${meta.run}_${gene}_avg_coverage_per_ccds_base_per_cell.tsv")

  script:
  """
  #!/usr/bin/env Rscript

  # libraries
  library(magrittr)

  # get gene info (incl. strand to orient 5'/3' ends)
  g <- 
    readr::read_tsv("${gene_bed}",
                    col_names = c("chr", "start", "end", "gene", "strand"))

  # create list of regions
  regions <- list()
  
  # get exonic regions
  regions[["exonic"]] <-
    readr::read_tsv("${exonic_bed}",
                    col_names = c("chr", "start", "end")) %>%
    # label each region, expand to all positions
    dplyr::distinct() %>%
    dplyr::mutate(arr = dplyr::case_when(g\$strand == "+" ~ start,
                                         g\$strand == "-" ~ -start)) %>%
    dplyr::arrange(arr) %>%
    dplyr::mutate(region = dplyr::row_number()) %>%
    dplyr::group_by(chr, region) %>%
    dplyr::reframe(pos = start:end) %>%
    # get exonic distance from 3'
    dplyr::group_by(pos) %>%
    dplyr::mutate(exonic_distance_from_5_prime = dplyr::cur_group_id()) %>%
    dplyr::ungroup()
  
  # get ccds regions
  regions[["ccds"]] <-
    readr::read_tsv("${ccds_bed}") %>%
    dplyr::mutate(arr = dplyr::case_when(g\$strand == "+" ~ start,
                                         g\$strand == "-" ~ -start)) %>%
    dplyr::arrange(arr) %>%
    dplyr::mutate(region = dplyr::row_number()) %>%
    dplyr::group_by(chr, region) %>%
    dplyr::reframe(pos = start:end)

  # get features (must wrangle from GFF3 attributes)
  features <-
    readr::read_tsv("${features_bed}") %>%
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

  # get mutations from Mitchell 2022, with calls in PB panel
  mut_annots <-
    readr::read_tsv("${mut_annots}") %>%
    dplyr::filter(gene == "${gene}")
  if (nrow(mut_annots) > 0) {
    mut_annots <-
      mut_annots %>%
      dplyr::right_join(
        tibble::tibble(chr = unique(features\$chr),
                       pos = min(features\$start):max(features\$end))) %>%
        dplyr::mutate(
          var_label = dplyr::case_when(alt_depth == 0 ~ NA,
                                      !is.na(HGVSp) ~ HGVSp,
                                      !is.na(HGVSc) ~ HGVSc,
                                      TRUE ~ mut_ref),
          label = dplyr::case_when(
            !is.na(var_label) ~ sprintf("%s (%s/%s)", var_label, alt_depth, depth),
            TRUE ~ NA))
  } else {
    mut_annots <- tibble::tibble()
  }

  # get coverage data
  cov <- readr::read_tsv("${cov}")

  # get all cell barcodes, remove suffix
  cell_barcodes <- readLines("${cell_barcodes}") %>% gsub("-1\$", "", .)

  # get celltype annotations
  celltypes <- readr::read_csv("${celltypes}")

  # get lists of celltype-annotated barcodes and list of all barcodes
  barcodes_ls <- split(celltypes\$barcode_${meta.kit}, celltypes\$celltype)
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

    }) %>%
    dplyr::bind_rows(.id = "celltype") %>%
    dplyr::mutate(celltype_label = paste0(celltype, " (", total_cells, ")") %>%
                    forcats::fct_reorder(total_cells, .desc = TRUE))

  # get total n cells (across celltypes)
  total_cells_all_cts <-
    cov_per_ct %>%
    dplyr::distinct(celltype_label, total_cells) %>%
    {sum(.\$total_cells)}

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
  
  # get average coverage per cell across the cCDS
  # TODO: add metric of coverage per base per expressing cell
  # expressing cells are cells with any coverage in the gene
  cov_per_ct %>%
    dplyr::right_join(regions[["ccds"]]) %>%
    dplyr::mutate(n_cells_x_cov = n_cells * coverage) %>%
    dplyr::summarise(total_cov = sum(n_cells_x_cov),
                     total_cells_x_bases = sum(n_cells),
                     avg_cov_per_base_per_cell = total_cov / total_cells_x_bases) %>%
    readr::write_tsv(
      "${meta.run}_${gene}_avg_coverage_per_ccds_base_per_cell.tsv")

  # save
  g %>% saveRDS("${meta.run}_${gene}.rds")
  regions %>% saveRDS("${meta.run}_${gene}_regions.rds")
  features %>% saveRDS("${meta.run}_${gene}_features.rds")
  mut_annots %>% saveRDS("${meta.run}_${gene}_mut_annots.rds")
  cov_per_ct %>% saveRDS("${meta.run}_${gene}_cov_per_ct.rds")
  exonic_cov %>%
    readr::write_tsv("${meta.run}_${gene}_coverage_per_exonic_position.tsv")
  """
}

// plot coverage
process plot_coverage {
  tag "${meta.run}_${gene}"
  label 'normal10gb'
  maxRetries 10
  publishDir "${params.out_dir}/runs/${meta.run}/${gene}/", mode: "copy"
  beforeScript "module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=~/R-tmp-4.4"

  input:
  tuple val(meta), val(gene), path(rdss)
  
  output:
  path("${meta.run}_${gene}_*_plot.${params.plot_device}")
  path("${meta.run}_${gene}_pie_mutations_plot.rds")

  script:
  """
  #!/usr/bin/env Rscript

  # libraries
  library(magrittr)
  library(ggplot2)
  library(patchwork)
  
  # read plot data
  g <- readRDS("${meta.run}_${gene}.rds")
  regions <- readRDS("${meta.run}_${gene}_regions.rds")
  features <- readRDS("${meta.run}_${gene}_features.rds")
  mut_annots <- readRDS("${meta.run}_${gene}_mut_annots.rds")
  cov_per_ct <- readRDS("${meta.run}_${gene}_cov_per_ct.rds")

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
        dplyr::select(dplyr::any_of(c("pos", "mut_type", "ref", "alt", "region"))) %>%
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
      geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
                aes(colour = coverage), position = "stack", width = 1,
                linewidth = 0.0001) +
      geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
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
      dplyr::summarise(`% of cells\\nw/ cov` = 100 *
                          sum(n_cells[coverage >= 2]) / sum(n_cells)) %>%
      ggplot(aes(x = pos, y = 1, fill = `% of cells\\nw/ cov`,
                  colour = `% of cells\\nw/ cov`)) +
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
      dplyr::filter(!is.na(mut_ref))
    mut_sites <-
      mut_sites %>%
      dplyr::group_by(mut_ref, chr, pos, gene) %>%
      dplyr::summarise(ref_depth = sum(ref_depth),
                      alt_depth = sum(alt_depth),
                      depth = sum(depth),
                      celltype = "all") %>%
      dplyr::bind_rows(mut_sites) %>%
      dplyr::mutate(
        region = mut_ref,
        radius = 0.45,
        celltype = factor(celltype, levels = ct_order),
        celltype_n = as.integer(celltype),
        arr = ifelse(g\$strand == "+", pos, -pos),
        mut_ref = forcats::fct_reorder(mut_ref, arr),
        mut_ref_n = as.integer(mut_ref),
        label = paste0(ifelse(depth == 0, "-", round(100 * alt_depth / depth, 1)),
                      "%\\n", alt_depth, "/", depth),
        `variant detected?` = ifelse(alt_depth > 0, "yes", NA)) %>%
      dplyr::select(-depth)

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
      scale_x_continuous(breaks = unique(mut_sites\$mut_ref_n),
                        labels = unique(mut_sites\$mut_ref),
                        guide = guide_axis(angle = -45),
                        expand = c(0, 0)) +
      scale_y_continuous(breaks = unique(mut_sites\$celltype_n),
                        labels = unique(mut_sites\$celltype),
                        expand = c(0, 0)) +
      theme_minimal() +
      scale_fill_manual(values = c("ref_depth" = "#c0d7fa",
                                  "alt_depth" = "#eabac7")) +
      scale_colour_manual(values = c("yes" = "red"), na.translate = FALSE, na.value = NA) +
      labs(x = "variant", title = p_title,
           subtitle = "Proportion of alt reads / all reads at the site") +
      theme(panel.grid = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_line())

    # plot ge, by celltype
    p_mut_ge <-
      cov_per_ct %>%
      dplyr::mutate(celltype = factor(celltype, levels = ct_order)) %>%
      dplyr::inner_join(mut_sites) %>%
      dplyr::group_by(mut_ref, celltype, `variant detected?`) %>%
      dplyr::summarise(n_cells_w_cov = sum(n_cells[coverage >= 1]),
                      n_cells_total = sum(n_cells),
                      `% of cells\\nw/ reads` = 100 * n_cells_w_cov / n_cells_total) %>%
      ggplot(aes(x = mut_ref, y = celltype)) +
      geom_tile(aes(fill = `% of cells\\nw/ reads`), width = 1) +
      geom_tile(aes(colour = `variant detected?`), fill = NA, linewidth = 1) +
      geom_label(aes(label = paste0(round(`% of cells\\nw/ reads`, 1), "%\\n",
                                    n_cells_w_cov, "/", n_cells_total))) +
      viridis::scale_fill_viridis() +
      scale_colour_manual(values = c("yes" = "red"), na.translate = FALSE, na.value = NA) +
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
    p_dist_polyA <-
      mut_sites %>%
      dplyr::mutate(pos_polyA = ifelse(g\$strand == "+", g\$end, g\$start),
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
    p_mut_pie / p_mut_ge / p_dist_polyA
  }

  # initialise list of lists plots
  p <- list()
  p[["genic"]] <- p[["exonic"]] <- p[["ccds"]] <- list()
  p_subtitle <- paste0("${meta.run} ${gene} (", g\$strand, " strand, ",
                       prettyNum(g\$end - g\$start, big.mark = ",",
                                 scientific = FALSE),
                       "bp, min coverage = ${params.min_cov})")

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
  if ("${meta.kit}" == "PB" & "${meta.seq_type}" == "panel" & nrow(mut_annots) > 0) {

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
    saveRDS(p_mut_pie, "${meta.run}_${gene}_pie_mutations_plot.rds")
    ggsave(
      "${meta.run}_${gene}_pie_mutations_plot.${params.plot_device}",
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
  if (g\$strand == "-") {
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
      paste0("${meta.run}_${gene}_", lvl, "_cov_plot.${params.plot_device}"),
      p1[["wgs_mutations"]] /
      p1[["mutations"]] /
      p1[["cov"]] /
      p1[["ge"]] /
      p1[["transcripts"]] +
        plot_layout(heights = c(p_wgs_muts_height, p_muts_height, 4, 0.5, 2)),
      dpi = 500)
    ggsave(
      paste0("${meta.run}_${gene}_", lvl, "_cov_pct_plot.${params.plot_device}"),
      p1[["wgs_mutations"]] /
      p1[["mutations"]] /
      p1[["cov_pct"]] /
      p1[["ge"]] /
      p1[["transcripts"]] +
        plot_layout(heights = c(p_wgs_muts_height, p_muts_height, 4, 0.5, 2)),
      dpi = 500)
  })
  """
}

// plot 5' drop-off
process plot_5_prime_dropoff {
  tag "${gene}"
  label 'normal10gb'
  maxRetries 10
  publishDir "${params.out_dir}/genes/${gene}/",
    mode: "copy",
    pattern: "*.{png,rds}"
  errorStrategy 'retry'
  beforeScript "module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=~/R-tmp-4.4"

  input:
  tuple val(gene),
        path(exonic_coverages)

  output:
  path("${gene}_exonic_coverage_plot.${params.plot_device}")
  path("${gene}_exonic_coverage_plot.rds")

  script:
  """
  #!/usr/bin/env Rscript

  # libraries
  library(magrittr)
  library(ggplot2)

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
    geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
             aes(colour = coverage), position = "stack", width = 1,
             linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
             aes(fill = coverage), position = "stack", width = 1) +
    viridis::scale_fill_viridis(na.value = "grey") +
    viridis::scale_colour_viridis(na.value = "grey") +
    facet_grid(run ~ ., scales = "free_y") +
    theme_minimal()
  
  # save
  ggsave("${gene}_exonic_coverage_plot.${params.plot_device}",
         p, dpi = 500)
  saveRDS(p, "${gene}_exonic_coverage_plot.rds")
  """
}

// main workflow
workflow {
  
  // get metadata + bam paths  
  Channel.fromPath(params.samplesheet, checkIfExists: true)
  | splitCsv(header: true)
  | map { row ->
      run = row.kit + "_" + row.seq_type + "_" + row.id
      meta = [id:row.id, seq_type:row.seq_type, kit:row.kit, run:run]
      [meta,
       file(row.bam, checkIfExists: true),
       file(row.celltypes, checkIfExists: true),
       file(row.mut_annots, checkIfExists: true),
       file(row.cell_barcodes, checkIfExists: true)]
  }
  | set { input }
  
  // download or locally link bams
  if (params.location == "irods") {
    // download input from irods
    ch_bam = irods(input)
  }
  else if (params.location == "local") {
    // input are locally available
    ch_bam = local(input)
  }

  // check bam is not truncated before proceeding
  ch_checked_bam = check_bam(ch_bam)
  
  // get gencode gff3
  gencode_gff3 = file(params.gencode_gff3)

  // get dndscv refcds
  refcds = file(params.refcds)

  // get driver gene coords
  ch_driver = Channel.fromPath(params.drivers).splitText().map{it -> it.trim()}
  get_driver_gene_ccds(ch_driver, refcds)
  get_driver_gene_coords(get_driver_gene_ccds.out, gencode_gff3)

  // output bed of all gene coords
  bed_hdr = Channel.of("chr\tstart\tend\tgene\tstrand\n")
  bed_hdr.concat(get_driver_gene_coords.out.gene_bed)
  | collectFile(name: params.out_dir + '/genes/genes.bed')

  // plot coverage per base per gene per cell
  ch_samples_x_drivers = \
    ch_checked_bam.combine(get_driver_gene_coords.out.gene_features)
  subset_bams_to_drivers_and_barcodes(ch_samples_x_drivers)
  get_coverage_per_cell(subset_bams_to_drivers_and_barcodes.out)
  wrangle_data(get_coverage_per_cell.out)
  plot_coverage(wrangle_data.out.plot_data)

  // plot 5' drop-off per gene
  ch_exonic_cov_per_gene = wrangle_data.out.exonic_coverage.groupTuple()
  plot_5_prime_dropoff(ch_exonic_cov_per_gene)
  
}
