#!/usr/bin/env nextflow

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
  tuple val(meta), val(bam), path(celltypes),
        path(mut_annots), path(wgs_mut_annots)
  
  output:
  tuple val(meta), path("${meta.run}.bam"), path("${meta.run}.bam.bai"),
        path(celltypes), path(mut_annots), path(wgs_mut_annots)
  
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
  tuple val(meta), val(bam), path(celltypes),
        path(mut_annots), path(wgs_mut_annots)

  output:
  tuple val(meta), path("${meta.run}.bam"), path("${meta.run}.bam.bai"),
        path(celltypes), path(mut_annots), path(wgs_mut_annots)

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
  tuple val(meta), path(bam), path(bai), path(celltypes),
        path(mut_annots), path(wgs_mut_annots)

  output:
  tuple val(meta), path(bam), path(bai), path(celltypes),
        path(mut_annots), path(wgs_mut_annots)

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

  input:
  tuple val(gene), path(ccds_bed)
  path(gencode_gff3)

  output:
  tuple val(gene),
        path(ccds_bed),
        path("${gene}.bed"),
        path("${gene}_features.bed"),
        path("${gene}_exonic_regions.bed")

  script:
  """
  # get the proteinc-coding gene coords from the gff3 file (chr start end gene strand)
  zcat ${gencode_gff3} \\
  | grep -w ${gene} \\
  | grep -w gene_type=protein_coding \\
  | awk -F"\\t" -v OFS="\\t" '\$3 == "gene" {print \$1, \$4, \$5, "${gene}", \$7}' \\
  > ${gene}.bed

  # get the transcript/exon coords from the gff3 file
  echo -e "chr\\tstart\\tend\\tgene\\ttype\\tattributes" > ${gene}_features.bed
  zcat ${gencode_gff3} \\
  | grep -w ${gene} \\
  | awk -F"\t" -v OFS="\t" '\$3 == "exon" || \$3 == "transcript" || \$3 == "start_codon" {print \$1, \$4, \$5, "${gene}", \$3, \$9}' \\
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

// subset the BAMs to only the driver regions
process subset_bams_to_drivers {
  tag "${meta.run}_${gene}"
  label 'normal'
  errorStrategy 'ignore'
  publishDir "${params.out_dir}/runs/${meta.run}/${gene}/",
    mode: "copy",
    pattern: "*_subset.{bam,bam.bai}"

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mut_annots), path(wgs_mut_annots),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)

  output:
  tuple val(meta),
        path("${meta.run}_subset.bam"), path("${meta.run}_subset.bam.bai"),
        path(celltypes),
        path(mut_annots), path(wgs_mut_annots),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)

  script:
  """
  module load samtools-1.19/python-3.12.0 
  samtools view -L $gene_bed -b $bam > ${meta.run}_subset.bam
  samtools index -@ ${task.cpus} ${meta.run}_subset.bam
  """
}

// split subsetted BAMs by cell, count coverage
process get_coverage_per_cell {
  tag "${meta.run}_${gene}"
  label 'week10gb'
  errorStrategy 'retry'
  publishDir "${params.out_dir}/runs/${meta.run}/${gene}/",
    mode: "copy",
    pattern: "*_coverage_per_cell.tsv"

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mut_annots), path(wgs_mut_annots),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)

  output:
  tuple val(meta), path("${meta.run}_${gene}_coverage_per_cell.tsv"),
        path(celltypes),
        path(mut_annots), path(wgs_mut_annots),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)

  script:
  """
  module load samtools-1.19/python-3.12.0 

  echo "getting cell barcodes"
  # get unique cell barcodes, trim CB:Z: prefix and -1 suffix
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

// plot coverage
process plot_coverage {
  tag "${meta.run}_${gene}"
  label 'normal10gb'
  maxRetries 10
  publishDir "${params.out_dir}/runs/${meta.run}/${gene}/",
    mode: "copy",
    pattern: "*{coverage_per_celltype.tsv,coverage_per_exonic_position.tsv,coverage_per_ccds_base_per_cell.tsv,.pdf,.png}"
  errorStrategy 'retry'
  beforeScript "module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=~/R-tmp-4.4"

  input:
  tuple val(meta), path(cov),
        path(celltypes),
        path(mut_annots), path(wgs_mut_annots),
        val(gene), path(ccds_bed), path(gene_bed), path(features_bed),
        path(exonic_bed)
  
  output:
  tuple val(gene),
        path("${meta.run}_${gene}_coverage_per_exonic_position.tsv"),
        emit: exonic_coverage
  path("${meta.run}_${gene}_coverage_per_celltype.tsv")
  path("${meta.run}_${gene}_*_plot.${params.plot_device}")

  script:
  def barcode_col = (meta.kit == "TX") ? "barcode" : (meta.kit == "PB") ? "barcode_revc" : ""
  """
  #!/usr/bin/env Rscript

  # libraries
  library(magrittr)
  library(ggplot2)
  library(patchwork)
  library(biomaRt)

  # get gene info (incl. strand to orient 5'/3' ends)
  gene <- 
    readr::read_tsv("${gene_bed}",
                    col_names = c("chr", "start", "end", "gene", "strand"))

  # get exonic regions
  exonic_regions <-
    readr::read_tsv("${exonic_bed}",
                    col_names = c("chr", "start", "end")) %>%
    # label each region, expand to all positions
    dplyr::distinct()
  # number facets by gene orientation
  if (gene\$strand == "+") {
    exonic_regions <- dplyr::arrange(exonic_regions, start)
  } else {
    exonic_regions <- dplyr::arrange(exonic_regions, desc(start))
  }
  exonic_regions <- 
    exonic_regions %>%
    dplyr::mutate(exonic_region = dplyr::row_number()) %>%
    dplyr::group_by(chr, exonic_region) %>%
    dplyr::reframe(pos = start:end)
  
  # get ccds regions
  ccds_regions <-
    readr::read_tsv("${gene}_ccds_regions.bed")
  if (gene\$strand == "+") {
    ccds_regions <- dplyr::arrange(ccds_regions, start)
  } else {
    ccds_regions <- dplyr::arrange(ccds_regions, desc(start))
  }
  ccds_regions <-
    ccds_regions %>%
    dplyr::mutate(ccds_region = dplyr::row_number()) %>%
    dplyr::group_by(chr, ccds_region) %>%
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
    dplyr::mutate(
      # move protein-coding to top of plot
      transcript_type = forcats::fct_relevel(transcript_type, "protein_coding"))

  # get celltype annotations
  celltypes <-
    readr::read_csv("${celltypes}") %>%
    dplyr::filter(barcode_${meta.kit} %in% colnames(cov))

  # get coverage data
  cov <- readr::read_tsv("${cov}")

  # get the counts of each unique value in each matrix row, per celltype
  cov_per_ct <-
    split(celltypes\$barcode_${meta.kit}, celltypes\$celltype) %>%
    purrr::map(function(CBs) {

      # get a matrix of reads in the celltype
      mat <- 
        cov %>%
        dplyr::select(dplyr::any_of(CBs)) %>%
        as.matrix()

      if (ncol(mat) > 0) {
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
      }

    }) %>%
    dplyr::bind_rows(.id = "celltype") %>%
    dplyr::mutate(celltype = paste0(celltype, " (", total_cells, ")") %>%
                    forcats::fct_reorder(total_cells, .desc = T))

  # save
  cov_per_ct %>%
    readr::write_tsv("${meta.run}_${gene}_coverage_per_celltype.tsv")

  # get total n cells (across celltypes)
  total_cells_all_cts <-
    cov_per_ct %>%
    dplyr::distinct(celltype, total_cells) %>%
    {sum(.\$total_cells)}

  # pile up coverage across celltypes at each exonic position
  exonic_cov <-
    cov_per_ct %>%
    dplyr::inner_join(exonic_regions) %>%
    # get number of cells per coverage per position
    dplyr::group_by(chr, pos, gene, exonic_region, coverage) %>%
    dplyr::summarise(n_cells = sum(n_cells)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_cells = total_cells_all_cts)

  # get exonic distance from the 3' end
  if (gene\$strand == "+") {
    exonic_cov <- dplyr::arrange(exonic_cov, dplyr::desc(pos))
  } else {
    exonic_cov <- dplyr::arrange(exonic_cov, pos)
  }
  exonic_cov <-
    exonic_cov %>%
    dplyr::group_by(pos) %>%
    dplyr::mutate(exonic_distance_from_3_prime = dplyr::cur_group_id()) %>%
    dplyr::ungroup()
  
  # save
  exonic_cov %>%
    readr::write_tsv("${meta.run}_${gene}_coverage_per_exonic_position.tsv")

  # get average coverage per cell across the cCDS
  cov_per_ct %>%
    dplyr::inner_join(ccds_regions) %>%
    dplyr::mutate(n_cells_x_cov = n_cells * coverage) %>%
    dplyr::summarise(total_cov = sum(n_cells_x_cov),
                     total_cells_x_bases = sum(n_cells),
                     avg_cov_per_base_per_cell = total_cov / total_cells_x_bases) %>%
    readr::write_tsv(
      "${meta.run}_${gene}_avg_coverage_per_ccds_base_per_cell.tsv")

  # initialise list of plots
  p <- list()
  p_theme <- theme_minimal()
  p_subtitle <- paste0("${meta.run} ${gene} (", gene\$strand, " strand, ",
                        prettyNum(gene\$end - gene\$start, big.mark = ",",
                                  scientific = FALSE),
                        "bp, min coverage = ${params.min_cov})")

  # plot wgs mutations
  wgs_mut_annots <-
    readr::read_tsv("${wgs_mut_annots}") %>%
    dplyr::filter(gene == "${gene}")
  
  if (nrow(wgs_mut_annots) > 0) {

      p[["wgs_mutations"]] <-
        wgs_mut_annots %>%
        ggplot(aes(x = pos, y = 0, colour = mut_type)) +
        geom_point() +
        scale_y_discrete(expand = c(0, 0)) +
        theme_void()
      p[["wgs_mutations_exonic"]] <-
        wgs_mut_annots %>%
        dplyr::right_join(exonic_regions) %>%
        ggplot(aes(x = pos, y = 1, colour = mut_type)) +
        geom_point(na.rm = TRUE) +
        facet_grid(. ~ exonic_region, scales = "free_x", space = "free_x") +
        theme_void() +
        scale_colour_discrete(na.translate = FALSE) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme(strip.text.x = element_blank(),
              panel.spacing = unit(0, "lines"),
              panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
              strip.background = element_rect(color = "grey", fill = NA,
                                              linewidth = 0.1, linetype = "solid"))
        
    } else {
      p[["wgs_mutations"]] <- ggplot() + theme_void()
      p[["wgs_mutations_exonic"]] <- ggplot() + theme_void()
    }

    # add mutation counts as subtitle
    p[["wgs_mutations"]] <-
      p[["wgs_mutations"]] +
      labs(subtitle = paste(nrow(mut_annots), "mutations called in WGS from Mitchell 2022"))
    exonic_vars <-
      mut_annots %>%
      dplyr::filter(Start_Position %in% exonic_regions\$pos)
    p[["wgs_mutations_exonic"]] <-
      p[["wgs_mutations_exonic"]] +
      labs(subtitle = paste(nrow(exonic_vars), "exonic mutations called in WGS from Mitchell 2022"))

  # plot mutations (if they exist)
  if ("${mut_annots}" != "NO_FILE") {
    mut_annots <-
      readr::read_tsv("${mut_annots}", comment = "#") %>%
      dplyr::filter(Hugo_Symbol == gene\$gene)
    
    if (nrow(mut_annots) > 0) {

      p[["mutations"]] <-
        mut_annots %>%
        ggplot(aes(x = Start_Position, y = 0, colour = Variant_Type)) +
        geom_point() +
        scale_y_discrete(expand = c(0, 0)) +
        theme_void()
      p[["mutations_exonic"]] <-
        mut_annots %>%
        dplyr::mutate(pos = Start_Position) %>%
        dplyr::right_join(exonic_regions) %>%
        #dplyr::arrange(is.na(Variant_Type)) %>%
        ggplot(aes(x = pos, y = 1, colour = Variant_Type)) +
        geom_point(na.rm = TRUE) +
        facet_grid(. ~ exonic_region, scales = "free_x", space = "free_x") +
        theme_void() +
        scale_colour_discrete(na.translate = FALSE) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme(strip.text.x = element_blank(),
              panel.spacing = unit(0, "lines"),
              panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
              strip.background = element_rect(color = "grey", fill = NA,
                                              linewidth = 0.1, linetype = "solid"))
        
    } else {
      p[["mutations"]] <- ggplot() + theme_void()
      p[["mutations_exonic"]] <- ggplot() + theme_void()
    }

    # add mutation counts as subtitle
    p[["mutations"]] <-
      p[["mutations"]] +
      labs(subtitle = paste(nrow(mut_annots), "mutations called"))
    exonic_vars <-
      mut_annots %>%
      dplyr::filter(Start_Position %in% exonic_regions\$pos)
    p[["mutations_exonic"]] <-
      p[["mutations_exonic"]] +
      labs(subtitle = paste(nrow(exonic_vars), "exonic mutations called"))
  
  } else {
    p[["mutations"]] <- ggplot() + theme_void()
    p[["mutations_exonic"]] <- ggplot() + theme_void()
  }

  # plot transcripts
  p[["transcripts"]] <-
    features %>%
    dplyr::mutate(exon_number = as.numeric(exon_number),
                  exon_alternating = as.character(exon_number %% 2)) %>%
    ggplot(aes(x = start, xend = end,
                y = transcript_id, yend = transcript_id, size = type)) +
    geom_segment(data = . %>% dplyr::filter(type != "start_codon"),
                  aes(colour = exon_alternating)) +
    scale_size_manual(values = c("transcript" = 0.5, "exon" = 3)) +
    scale_colour_manual(values = c("1" = "royalblue4", "0" = "steelblue3")) +
    # plot the start codon
    geom_point(data = . %>% dplyr::filter(type == "start_codon"),
                size = 3, colour = "lightskyblue1") +
    geom_point(data = . %>% dplyr::filter(type == "start_codon"),
                shape = 83, size = 2, colour = "blue4") +
    facet_grid(transcript_type ~ ., space = "free_y", scales = "free_y") +
    p_theme +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          axis.text.x = element_text(angle = -90, hjust = 1),
          axis.text.y = element_blank()) +
    guides(colour = "none") +
    labs(x = "position")

  # plot transcripts (exonic only)
  p[["transcripts_exonic"]] <-
    features %>%
    dplyr::filter(type == "exon") %>%
    dplyr::mutate(exon_number = as.numeric(exon_number),
                  exon_alternating = as.character(exon_number %% 2)) %>%
    # get exonic facets
    dplyr::group_by(dplyr::across(c(-start, -end))) %>%
    dplyr::reframe(pos = start:end) %>%
    dplyr::inner_join(exonic_regions) %>%
    dplyr::group_by(dplyr::across(c(-pos))) %>%
    dplyr::summarise(start = min(pos), end = max(pos)) %>%
    # plot
    ggplot(aes(x = start, xend = end,
                y = transcript_id, yend = transcript_id, size = type)) +
    geom_segment(data = . %>% dplyr::filter(type != "start_codon"),
                  aes(colour = exon_alternating)) +
    scale_size_manual(values = c("transcript" = 0.5, "exon" = 3)) +
    scale_colour_manual(values = c("1" = "royalblue4", "0" = "steelblue3")) +
    scale_x_discrete(expand = c(0, 0)) +
    # plot the start codon
    geom_point(data = . %>% dplyr::filter(type == "start_codon"),
                size = 3, colour = "lightskyblue1") +
    geom_point(data = . %>% dplyr::filter(type == "start_codon"),
                shape = 83, size = 2, colour = "blue4") +
    facet_grid(transcript_type ~ exonic_region, space = "free", scales = "free") +
    p_theme +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          strip.text.x = element_blank(),
          axis.text.x = element_text(angle = -90, hjust = 1, size = 2),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
          strip.background = element_rect(color = "grey", fill = NA,
                                          linewidth = 0.1, linetype = "solid")) +
    guides(colour = "none") +
    labs(x = "position")
  
  # plot transcripts (ccds only)
  p[["transcripts_ccds"]] <-
    features %>%
    dplyr::filter(type == "exon") %>%
    dplyr::mutate(exon_number = as.numeric(exon_number),
                  exon_alternating = as.character(exon_number %% 2)) %>%
    # get exonic facets
    dplyr::group_by(dplyr::across(c(-start, -end))) %>%
    dplyr::reframe(pos = start:end) %>%
    dplyr::inner_join(ccds_regions) %>%
    dplyr::group_by(dplyr::across(c(-pos))) %>%
    dplyr::summarise(start = min(pos), end = max(pos)) %>%
    # plot
    ggplot(aes(x = start, xend = end,
                y = transcript_id, yend = transcript_id, size = type)) +
    geom_segment(data = . %>% dplyr::filter(type != "start_codon"),
                  aes(colour = exon_alternating)) +
    scale_size_manual(values = c("transcript" = 0.5, "exon" = 3)) +
    scale_colour_manual(values = c("1" = "royalblue4", "0" = "steelblue3")) +
    scale_x_discrete(expand = c(0, 0)) +
    # plot the start codon
    geom_point(data = . %>% dplyr::filter(type == "start_codon"),
                size = 3, colour = "lightskyblue1") +
    geom_point(data = . %>% dplyr::filter(type == "start_codon"),
                shape = 83, size = 2, colour = "blue4") +
    facet_grid(transcript_type ~ ccds_region, space = "free", scales = "free") +
    p_theme +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          strip.text.x = element_blank(),
          axis.text.x = element_text(angle = -90, hjust = 1, size = 2),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
          strip.background = element_rect(color = "grey", fill = NA,
                                          linewidth = 0.1, linetype = "solid")) +
    guides(colour = "none") +
    labs(x = "position")

  # plot coverage per base
  p[["cov"]] <-
    cov_per_ct %>%
    dplyr::arrange(-coverage) %>%
    ggplot(aes(x = pos, y = n_cells)) +
    # plot grey background (0 coverage)
    geom_col(fill = NA, colour = "grey", position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(fill = "grey86", position = "stack", width = 1) +
    # plot dark grey background (0 < cov < min_cov)
    geom_col(data = . %>% dplyr::filter(coverage > 0, coverage < ${params.min_cov}),
              fill = NA, colour = "darkgrey", position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage > 0, coverage < ${params.min_cov}),
              fill = "grey48", position = "stack", width = 1) +
    # plot coverage
    geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
              aes(colour = coverage), position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
              aes(fill = coverage), position = "stack", width = 1) +
    viridis::scale_fill_viridis(na.value = "grey") +
    facet_grid(celltype ~ ., scales = "free_y") +
    p_theme +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 3),
          axis.title.x = element_blank()) +
    labs(title = "Coverage per cell", subtitle = p_subtitle, y = "n cells") +
    guides(colour = "none")

  # plot coverage per base (exonic only)
  p[["cov_exonic"]] <-
    cov_per_ct %>%
    dplyr::inner_join(exonic_regions) %>%
    dplyr::arrange(-coverage) %>%
    ggplot(aes(x = pos, y = n_cells)) +
    # plot grey background (0 coverage)
    geom_col(fill = NA, colour = "grey", position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(fill = "grey86", position = "stack", width = 1) +
    # plot dark grey background (0 < cov < min_cov)
    geom_col(data = . %>% dplyr::filter(coverage > 0, coverage < ${params.min_cov}),
              fill = NA, colour = "darkgrey", position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage > 0, coverage < ${params.min_cov}),
              fill = "grey48", position = "stack", width = 1) +
    # plot coverage
    geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
              aes(colour = coverage), position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
              aes(fill = coverage), position = "stack", width = 1) +
    viridis::scale_fill_viridis(na.value = "grey") +
    scale_x_discrete(expand = c(0, 0)) +
    facet_grid(celltype ~ exonic_region, scales = "free", space = "free_x") +
    p_theme +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          strip.text.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 3),
          axis.title.x = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
          strip.background = element_rect(color = "grey", fill = NA,
                                          linewidth = 0.1, linetype = "solid")) +
    labs(title = "Exonic coverage per cell", subtitle = p_subtitle, y = "n cells") +
    guides(colour = "none")

  # plot coverage per base (exonic only)
  p[["cov_ccds"]] <-
    cov_per_ct %>%
    dplyr::inner_join(ccds_regions) %>%
    dplyr::arrange(-coverage) %>%
    ggplot(aes(x = pos, y = n_cells)) +
    # plot grey background (0 coverage)
    geom_col(fill = NA, colour = "grey", position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(fill = "grey86", position = "stack", width = 1) +
    # plot dark grey background (0 < cov < min_cov)
    geom_col(data = . %>% dplyr::filter(coverage > 0, coverage < ${params.min_cov}),
              fill = NA, colour = "darkgrey", position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage > 0, coverage < ${params.min_cov}),
              fill = "grey48", position = "stack", width = 1) +
    # plot coverage
    geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
              aes(colour = coverage), position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage >= ${params.min_cov}),
              aes(fill = coverage), position = "stack", width = 1) +
    viridis::scale_fill_viridis(na.value = "grey") +
    scale_x_discrete(expand = c(0, 0)) +
    facet_grid(celltype ~ ccds_region, scales = "free", space = "free_x") +
    p_theme +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          strip.text.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 3),
          axis.title.x = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
          strip.background = element_rect(color = "grey", fill = NA,
                                          linewidth = 0.1, linetype = "solid")) +
    labs(title = "Canonical CDS coverage per cell", subtitle = p_subtitle, y = "n cells") +
    guides(colour = "none")

  # plot overall genotyping efficiency below celltype breakdown
  p[["ge"]] <-
    cov_per_ct %>%
    dplyr::group_by(pos) %>%
    dplyr::summarise(`% of cells\\nw/ cov` = 100 *
                        sum(n_cells[coverage >= ${params.min_cov}]) / sum(n_cells)) %>%
    ggplot(aes(x = pos, y = 1, fill = `% of cells\\nw/ cov`,
                colour = `% of cells\\nw/ cov`)) +
    geom_tile(width = 1, linewidth = 0.001) +
    viridis::scale_fill_viridis() +
    viridis::scale_colour_viridis() +
    p_theme +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  # plot overall genotyping efficiency below celltype breakdown (exonic only)
  p[["ge_exonic"]] <-
    cov_per_ct %>%
    dplyr::inner_join(exonic_regions) %>%
    dplyr::group_by(exonic_region, pos) %>%
    dplyr::summarise(`% of cells\\nw/ cov` = 100 *
                        sum(n_cells[coverage >= ${params.min_cov}]) / sum(n_cells)) %>%
    ggplot(aes(x = pos, y = 1, fill = `% of cells\\nw/ cov`,
                colour = `% of cells\\nw/ cov`)) +
    geom_tile(width = 1, linewidth = 0.001) +
    viridis::scale_fill_viridis() +
    viridis::scale_colour_viridis() +
    scale_x_discrete(expand = c(0, 0)) +
    facet_grid(. ~ exonic_region, scales = "free", space = "free_x") +
    p_theme +
    theme(strip.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
          strip.background = element_rect(color = "grey", fill = NA,
                                          linewidth = 0.1, linetype = "solid"))

  # plot genotyping efficiency (% of cells with >= `min_cov` reads per position)
  p[["ge_per_ct"]] <-
    cov_per_ct %>%
    dplyr::group_by(celltype, pos, total_cells) %>%
    dplyr::summarise(`% of cells\\nw/ cov` =
                        100 *
                        sum(n_cells[coverage >= ${params.min_cov}]) /
                        sum(n_cells)) %>%
    ggplot(aes(x = pos, y = total_cells, fill = `% of cells\\nw/ cov`,
                colour = `% of cells\\nw/ cov`)) +
    geom_col(width = 1, linewidth = 0.005) +
    viridis::scale_fill_viridis() +
    viridis::scale_colour_viridis() +
    facet_grid(celltype ~ ., space = "free_y", scales = "free_y") +
    p_theme +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off",
          axis.text.y = element_blank(), axis.title.x = element_blank()) +
    labs(title = "Genotyping efficiency per celltype",
          subtitle = p_subtitle)

  # plot average coverage per base
  p[["avg_cov"]] <-
    cov_per_ct %>%
    dplyr::group_by(celltype, pos) %>%
    dplyr::summarise(mean_coverage = mean(coverage),
                      median_coverage = median(coverage)) %>%
    ggplot(aes(x = pos, y = median_coverage)) +
    geom_line() +
    facet_grid(celltype ~ ., scales = "free_y") +
    p_theme +
    theme(strip.text.y.right = element_text(angle = 0), strip.clip = "off") +
    labs(title = "Median coverage per celltype per position",
          subtitle = p_subtitle,
          x = "position", y = "coverage")

  # if strand is "-", flip the x axis
  if (gene\$strand == "-") {
    p <- lapply(p, function(x) x + scale_x_reverse(expand = c(0, 0)))
  }

  # save plots
  ggsave("${meta.run}_${gene}_coverage_plot.${params.plot_device}",
        p[["wgs_mutations"]] / p[["mutations"]] / p[["cov"]] / p[["ge"]] /
        p[["transcripts"]] +
          plot_layout(heights = c(0.5, 0.5, 4, 0.5, 2)),
        dpi = 500)

  ggsave("${meta.run}_${gene}_coverage_exonic_plot.${params.plot_device}",
        p[["wgs_mutations_exonic"]] / p[["mutations_exonic"]] /
        p[["cov_exonic"]] / p[["ge_exonic"]] / p[["transcripts_exonic"]] +
          plot_layout(heights = c(1, 4, 0.5, 2)),
        dpi = 500)

  ggsave("${meta.run}_${gene}_genotyping_efficiency_plot.${params.plot_device}",
        p[["ge_per_ct"]] / p[["transcripts"]] + plot_layout(heights = c(4, 1)),
        dpi = 500)

  ggsave("${meta.run}_${gene}_coverage_median_plot.${params.plot_device}",
        p[["avg_cov"]] / p[["transcripts"]] + plot_layout(heights = c(4, 1)),
          dpi = 500)
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
    ggplot(aes(x = exonic_distance_from_3_prime, y = n_cells)) +
    # plot grey background (0 coverage)
    geom_col(fill = NA, colour = "grey", position = "stack", width = 1,
             linewidth = 0.0001) +
    geom_col(fill = "grey86", position = "stack", width = 1) +
    # plot dark grey background (0 < cov < min_cov) +
    geom_col(data = . %>% dplyr::filter(coverage > 0, coverage < ${params.min_cov}),
              fill = NA, colour = "darkgrey", position = "stack", width = 1,
              linewidth = 0.0001) +
    geom_col(data = . %>% dplyr::filter(coverage > 0, coverage < ${params.min_cov}),
              fill = "grey48", position = "stack", width = 1) +
    # plot coverage
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
      if (row.mut_annots == "NA") {
        mut_annots = file("${projectDir}/assets/NO_FILE", checkIfExists: true) 
      } else {
        mut_annots = file(row.mut_annots, checkIfExists: true)
      }
      [meta,
       file(row.bam, checkIfExists: true),
       file(row.celltypes, checkIfExists: true),
       mut_annots,
       file(row.wgs_mut_annots, checkIfExists: true)]
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

  // plot coverage per base per gene per cell
  ch_samples_x_drivers = ch_checked_bam.combine(get_driver_gene_coords.out)
  subset_bams_to_drivers(ch_samples_x_drivers)
  get_coverage_per_cell(subset_bams_to_drivers.out)
  plot_coverage(get_coverage_per_cell.out)

  // plot 5' drop-off per gene
  ch_exonic_cov_per_gene = plot_coverage.out.exonic_coverage.groupTuple()
  plot_5_prime_dropoff(ch_exonic_cov_per_gene)
  
}
