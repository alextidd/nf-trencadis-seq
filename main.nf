#!/usr/bin/env nextflow
// TODO: simplify mutations + fix handling in plot_mut_pie()
// TODO: fold Rmd reports into the nextflow pipeline
// TODO: add gene expression distribution alongside mutations / coverage plots
// TODO: create summary plot per gene
// TODO: add celltype marker genes to the expression dotplot (ask Laura)
// TODO: check over the transcript IDs that Carol from pipeline sent - begin deconvolution in the dataset

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { get_irods_bam  } from '../modules/get_irods_bam.nf'
include { get_local_bam  } from '../modules/get_local_bam.nf'
include { samtools_index } from '../modules/samtools_index.nf'
include { MOSDEPTH       } from './modules/nf-core/mosdepth/main'
include { single_cell    } from './modules/local/single_cell.nf'
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// check bam is not truncated before proceeding
process check_bam {
  tag "${meta.id}"
  label 'normal'

  input:
  tuple val(meta), path(bam),
        path(mutations), path(celltypes), path(cell_barcodes)

  output:
  tuple val(meta),
        path(bam), path("${meta.id}.${meta.bam_ext}.${meta.bai_ext}"),
        path(mutations), path(celltypes), path(cell_barcodes)

  script:
  """
  # modules
  module load samtools-1.19/python-3.12.0 
  samtools quickcheck ${bam}

  # index
  samtools index -@ ${task.cpus} ${bam}
  """
}

// get gene coords
process get_gene_coords {
  tag "${gene}"
  label 'normal10gb'
  publishDir "${params.out_dir}/genes/${gene}/",
    mode: "copy",
    pattern: "*.bed"

  input:
  val(gene)
  path(gencode_gff3)
  path(refcds)

  output:
  tuple val(gene),
        path("${gene}.bed"),
        path("${gene}_regions.rds"),
        path("${gene}_features.rds"), emit: gene_features
  path "${gene}.bed", emit: gene_bed

  script:
  def no_chr = params.no_chr ? "TRUE" : "FALSE"
  """
  get_gene_coords.R \\
    --gene ${gene} \\
    --gencode_gff3 ${gencode_gff3} \\
    --refcds ${refcds} \\
    --no_chr ${no_chr}
  """
}

// collate gene coords
process collate_gene_coords {
  label 'normal'
  publishDir "${params.out_dir}/genes/", mode: "copy"

  input:
  path(genes)

  output:
  path("genes.bed")

  script:
  """
  echo -e "chr\tstart\tend\tgene\tstrand" > genes.bed
  cat ${genes} >> genes.bed
  """
}

// subset the BAM to only reads in the genes
process subset_bam_to_genes {
  tag "${meta.id}_${meta.gene}"
  label 'normal4core'

  input:
  tuple val(meta), path(bam), path(bai),
        path(mutations), path(celltypes), path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features)

  output:
  tuple val(meta),
        path("${meta.id}_subset_genes.bam"),
        path("${meta.id}_subset_genes.bam.bai"),
        path(mutations), path(celltypes), path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features)

  script:
  """
  module load samtools-1.19/python-3.12.0 

  # subset to reads that fall in gene regions
  samtools view -@ ${task.cpus} -L ${gene_bed} -b ${bam} \\
  > ${meta.id}_subset_genes.bam

  # index
  samtools index -@ ${task.cpus} ${meta.id}_subset_genes.bam
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
  
  // validate input parameters
  validateParameters()

  // print summary of supplied parameters
  log.info paramsSummaryLog(workflow)

  // get input bams
  Channel
    .fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
  | map { meta, bam, mutations, celltypes, cell_barcodes ->
          [meta, bam]
  }
    | set { in_bam }

  // get samplesheet
  Channel
  .fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
  | map { meta, bam, mutations, celltypes, cell_barcodes ->
          [meta, mutations, celltypes, cell_barcodes]
  }
  | set { input }
  
  // download or locally link bams
  if (params.location == "irods") {
    // download input from irods
    ch_bam = get_irods_bam(in_bam)
  } else if (params.location == "local") {
    // inputs are locally available
    ch_bam = get_local_bam(in_bam)
  }

  // index bams
  samtools_index(ch_bam)
  
  // get gencode gff3
  gencode_gff3 = file(params.gencode_gff3, checkIfExists: true)

  // get refcds
  refcds = file(params.refcds, checkIfExists: true)

  // get gene coords
  ch_gene = Channel.fromPath(params.genes).splitText().map{it -> it.trim()}
  get_gene_coords(ch_gene, gencode_gff3, refcds)

  // output bed of all gene coords
  get_gene_coords.out.gene_bed.collect()
  | collate_gene_coords

  // combine gene-x-id, subset bams
  samtools_index.out \
  .join(input) \
  .combine(get_gene_coords.out.gene_features)
  | map { meta, bam, bai, mutations, celltypes, cell_barcodes,
          gene, gene_bed, gene_features, exonic_bed ->
          [meta + [gene: gene],
            bam, bai, mutations, celltypes, cell_barcodes, 
            gene_bed, gene_features, exonic_bed]
  }
  | subset_bam_to_genes
  | set { ch_bam_x_gene }

  // get coverage per celltype
  single_cell(ch_bam_x_gene)
  ch_cov_per_ct = single_cell.out.cov_per_ct
  ch_cell_bams = single_cell.out.cell_bams
  
  // // plot 5' drop-off per gene
  // ch_exonic_cov_per_gene = get_coverage_per_celltype.out.exonic_coverage.groupTuple()
  // plot_5_prime_dropoff(ch_exonic_cov_per_gene)
  
}