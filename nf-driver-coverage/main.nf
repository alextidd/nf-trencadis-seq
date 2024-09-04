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
  tuple val(meta), val(bam)
  
  output:
  tuple val(meta), path("${meta.run}.bam"), path("${meta.run}.bam.bai"), emit: bams
  
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
  tuple val(meta), val(bam)

  output:
  tuple val(meta), path("${meta.run}.bam"), path("${meta.run}.bam.bai")

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
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path(bam), path(bai)

  script:
  """
  module load samtools-1.19/python-3.12.0 
  samtools quickcheck $bam
  """
}

// get driver gene coords
process get_driver_gene_coords {
  tag "${gene}"
  label 'normal'

  input:
  val(gene)

  output:
  tuple val(gene), path("${gene}.bed"), path("${gene}_exons.bed")

  script:
  """
  # get the gene coords from the gff3 file (chr start end gene strand)
  zcat /lustre/scratch125/casm/team268im/at31/reference/gencode/gencode.v28.annotation.gff3.gz \
  | grep -w ${gene} \
  | awk -F"\\t" -v OFS="\\t" '\$3 == "gene" {print \$1, \$4, \$5, "${gene}", \$7}' \
  > ${gene}.bed

  # get the transcript/exon coords from the gff3 file
  echo -e "chr\\tstart\\tend\\tgene\\ttype\\tattributes" > ${gene}_exons.bed
  zcat /lustre/scratch125/casm/team268im/at31/reference/gencode/gencode.v28.annotation.gff3.gz \
  | grep -w ${gene} \
  | awk -F"\t" -v OFS="\t" '\$3 == "exon" || \$3 == "transcript" || \$3 == "start_codon" {print \$1, \$4, \$5, "${gene}", \$3, \$9}' \
  >> ${gene}_exons.bed
  """
}

// subset the BAMs to only the driver regions
process subset_bams_to_drivers {
  tag "${meta.run}_${gene}"
  label 'normal'
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(bam), path(bai),
        val(gene), path(gene_bed), path(exons_bed)

  output:
  tuple val(meta),
        path("${meta.run}_subset.bam"), path("${meta.run}_subset.bam.bai"),
        val(gene), path(gene_bed), path(exons_bed)

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
  label 'normal'
  errorStrategy 'ignore'
  publishDir "${params.out_dir}/${meta.run}/${gene}/",
    mode: "copy",
    pattern: "*_coverage_per_cell.tsv"

  input:
  tuple val(meta), path(bam), path(bai),
        val(gene), path(gene_bed), path(exons_bed)

  output:
  tuple val(meta), path("${meta.run}_${gene}_coverage_per_cell.tsv"),
        val(gene), path(gene_bed), path(exons_bed)

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
    mkdir \${gene}
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

  echo "cleaning up cell_bams/*"
  rm -rf cell_bams

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
  publishDir "${params.out_dir}/${meta.run}/${gene}/",
    mode: "copy",
    pattern: "*{coverage_per_celltype.tsv,.pdf,.png}"
  errorStrategy 'retry'
  beforeScript "module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=~/R-tmp-4.4"

  input:
  tuple val(meta), path(cov), val(gene), path(gene_bed), path(exons_bed)
  path(celltypes)
  
  output:
  path("*_plot.${params.plot_device}")
  path("${meta.run}_${gene}_coverage_per_celltype.tsv")

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

  # get canonical transcripts for genes
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  canonical_transcript <-
    getBM(attributes = c("ensembl_transcript_id", "transcript_is_canonical"),
          filters = "hgnc_symbol",
          values = "${gene}",
          mart = mart) %>%
    dplyr::filter(transcript_is_canonical == 1) %>%
    dplyr::pull(ensembl_transcript_id)

  # get exons (must wrangle from GFF3 attributes)
  exons <-
    readr::read_tsv("${exons_bed}") %>%
    dplyr::mutate(
      id = gsub(";.*", "", attributes)) %>%
    tidyr::separate_longer_delim(cols = "attributes", delim = ";") %>%
    tidyr::separate_wider_delim(
      cols = "attributes", delim = "=", names = c("name", "value")) %>%
    dplyr::filter(name %in% c("exon_id", "transcript_id", "transcript_type",
                              "exon_number")) %>%
    tidyr::pivot_wider() %>%
    dplyr::mutate(
      transcript_id = gsub("\\\\..*", "", transcript_id),
      # mark canonical, move canonical + protein-coding to top of plot
      transcript_type = dplyr::case_when(
        transcript_id == canonical_transcript ~ paste0("canonical_",
                                                       transcript_type),
        TRUE ~ transcript_type) %>%
      forcats::fct_relevel("canonical_protein_coding", "protein_coding"))

  # get coverage data, remove `-1` suffix from colnames
  cov <-
    readr::read_tsv("${cov}")
  colnames(cov) <- colnames(cov) %>% stringr::str_remove("-1\$")

  # get celltype annotations
  annot <- readr::read_tsv("${celltypes}") %>%
    dplyr::filter(id == "${meta.id}") %>%
    dplyr::select(CB = ${barcode_col}, celltype) %>%
    dplyr::filter(CB %in% colnames(cov))

  # get the counts of each unique value in each matrix row, per celltype
  cov_per_ct <-
    split(annot\$CB, annot\$celltype) %>%
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

  # initialise list of plots
  p <- list()
  p_theme <- theme_minimal()
  p_subtitle <- paste0("${meta.run} ${gene} (", gene\$strand, " strand, ",
                        prettyNum(gene\$end - gene\$start, big.mark = ",",
                                  scientific = FALSE),
                        "bp, min coverage = ${params.min_cov})")

  # plot protein-coding transcripts
  p[["transcripts"]] <-
    exons %>%
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
          axis.text.y = element_blank()) +
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
          axis.text.y = element_text(size = 3),
          axis.title.x = element_blank()) +
    labs(title = "Coverage per cell", subtitle = p_subtitle, y = "n cells") +
    guides(colour = "none")
  
  # plot overall genotyping efficiency below celltype breakdown
  p[["ge"]] <-
    cov_per_ct %>%
    dplyr::group_by(pos) %>%
    dplyr::summarise(`genotyping\\nefficiency (%)` = 100 *
                        sum(n_cells[coverage >= 2]) / sum(n_cells)) %>%
    ggplot(aes(x = pos, y = 1, fill = `genotyping\\nefficiency (%)`,
                colour = `genotyping\\nefficiency (%)`)) +
    geom_tile(width = 1, linewidth = 0.001) +
    viridis::scale_fill_viridis() +
    viridis::scale_colour_viridis() +
    p_theme +
    theme(axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  # plot genotyping efficiency (% of cells with >= `min_cov` reads per position)
  p[["ge_per_ct"]] <-
    cov_per_ct %>%
    dplyr::group_by(celltype, pos, total_cells) %>%
    dplyr::summarise(`genotyping\\nefficiency (%)` =
                        100 *
                        sum(n_cells[coverage >= ${params.min_cov}]) /
                        sum(n_cells)) %>%
    ggplot(aes(x = pos, y = total_cells, fill = `genotyping\\nefficiency (%)`,
                colour = `genotyping\\nefficiency (%)`)) +
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
    p <- lapply(p, function(x) x + scale_x_reverse())
  }

  # save plots
  ggsave("${meta.run}_${gene}_coverage_plot.${params.plot_device}",
         p[["cov"]] / p[["ge"]] / p[["transcripts"]] + plot_layout(heights = c(4, 0.5, 2)),
         dpi = 200)
  
  ggsave("${meta.run}_${gene}_genotyping_efficiency_plot.${params.plot_device}",
         p[["ge_per_ct"]] / p[["transcripts"]] + plot_layout(heights = c(4, 1)),
         dpi = 200)

  ggsave("${meta.run}_${gene}_median_coverage_plot.${params.plot_device}",
         p[["avg_cov"]] / p[["transcripts"]] + plot_layout(heights = c(4, 1)),
          dpi = 200)
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
      [meta, file(row.bam, checkIfExists: true)]
  }
  | set { samplesheet }
  
  // download or locally link bams
  if (params.location == "irods") {
    // download samplesheet from irods
    ch_bam = irods(samplesheet)
  }
  else if (params.location == "local") {
    // samplesheet are locally available
    ch_bam = local(samplesheet)
  }

  // check bam is not truncated before proceeding
  ch_checked_bam = check_bam(ch_bam)

  // get celltypes
  celltypes = file(params.celltypes)
  
  // get driver gene coords
  ch_driver = Channel.fromPath(params.drivers).splitText().map{it -> it.trim()}
  get_driver_gene_coords(ch_driver)

  // plot coverage per base per gene per cell
  ch_samples_x_drivers = ch_checked_bam.combine(get_driver_gene_coords.out)
  subset_bams_to_drivers(ch_samples_x_drivers)
  get_coverage_per_cell(subset_bams_to_drivers.out)
  plot_coverage(get_coverage_per_cell.out,
                celltypes)
  
}
