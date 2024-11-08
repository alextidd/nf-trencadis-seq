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
  wrangle_data.R \\
    --gene ${gene} \\
    --gene_bed ${gene_bed} \\
    --exonic_bed ${exonic_bed} \\
    --ccds_bed ${ccds_bed} \\
    --features_bed ${features_bed} \\
    --mut_annots ${mut_annots} \\
    --cov ${cov} \\
    --cell_barcodes ${cell_barcodes} \\
    --celltypes ${celltypes} \\
    --meta_kit ${meta.kit} \\
    --meta_run ${meta.run}
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
  path("${meta.run}_${gene}_pie_mutations_plot.rds"), optional: true

  script:
  """
  plot_coverage.R \\
    --meta_run ${meta.run} \\
    --meta_seq_type ${meta.seq_type} \\
    --meta_kit ${meta.kit} \\
    --gene ${gene} \\
    --min_cov ${params.min_cov} \\
    --plot_device ${params.plot_device}
  """
}

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
  plot_5_prime_dropoff.R \\
    --gene ${gene} \\
    --min_cov ${params.min_cov} \\
    --plot_device ${params.plot_device}
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
