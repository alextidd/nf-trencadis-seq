#!/usr/bin/env nextflow
// TODO: add genotyping step to pipeline
// TODO: simplify mutations + fix handling in plot_mut_pie()
// TODO: ignore noncoding mutations!
// TODO: fold Rmd reports into the nextflow pipeline
// TODO: fix called mutation annotations to fit + facet by celltype
// TODO: add gene expression distribution alongside mutations / coverage plots
// TODO: create summary plot per gene
// TODO: add celltype marker genes to the expression dotplot (ask Laura)
// TODO: check over the transcript IDs that Carol from pipeline sent - begin deconvolution in the dataset

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { single_cell } from './modules/single_cell.nf'
include { bulk        } from './modules/bulk.nf'

// Download a given sample's BAM from iRODS
// Then either retrieve the BAI or make one via indexing
// The maxForks of 10 was set after asking jc18 about best iRODS practices
process irods {
  tag "${meta.id}"
  maxForks 10
  label 'normal4core'

  input:
  tuple val(meta), val(bam),
        path(mutations), path(celltypes), path(cell_barcodes)
  
  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"),
        path(mutations), path(celltypes), path(cell_barcodes)
  
  script:
  """
  iget -K ${bam} ${meta.id}.bam
  if [[ `ils ${bam}.bai | wc -l` == 1 ]]
  then
      iget -K ${bam}.bai ${meta.id}.bam.bai
  else
      samtools index -@ ${task.cpus} ${meta.id}.bam
  fi
  """
}

// The equivalent of an irods download, but for a local copy of samplesheet
// Symlink the BAM/BAI appropriately so they're named the right thing for downstream
process local {
  tag "${meta.id}"
  maxForks 10
  label 'normal4core'

  input:
  tuple val(meta), val(bam),
        path(mutations), path(celltypes), path(cell_barcodes)

  output:
  tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"),
        path(mutations), path(celltypes), path(cell_barcodes)

  script:
  """
  # create local symbolic link 
  ln -s ${bam} ${meta.id}.bam
  if [ -f "${bam}.bai" ] ; then
      ln -s ${bam}.bai ${meta.id}.bam.bai
  else
      samtools index -@ ${task.cpus} ${meta.id}.bam
  fi
  """
}

// check bam is not truncated before proceeding
process check_bam {
  tag "${meta.id}"
  label 'normal'
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(bam), path(bai),
        path(mutations), path(celltypes), path(cell_barcodes)

  output:
  tuple val(meta), path(bam), path(bai),
        path(mutations), path(celltypes), path(cell_barcodes)

  script:
  """
  module load samtools-1.19/python-3.12.0 
  samtools quickcheck ${bam}
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
  """
  get_gene_coords.R \\
    --gene ${gene} \\
    --gencode_gff3 ${gencode_gff3} \\
    --refcds ${refcds}
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
  label 'normal'

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
  samtools view -L ${gene_bed} -b ${bam} > ${meta.id}_subset_genes.bam

  # index
  samtools index -@ ${task.cpus} ${meta.id}_subset_genes.bam
  """
}

// genotype the mutations
process genotype_mutations {
  tag "${meta.id}_${meta.gene}"
  label 'normal10gb'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy"
  
  input:
  tuple val(meta),
        path(mutations), path(celltypes),
        path(gene_bed),
        path(cell_bams, stageAs: "cell_bams/*")
  
  output:
  tuple val(meta),
        path("genotyped_mutations_per_celltype.tsv"), emit: geno_per_ct
  path("genotyped_mutations_per_cell.tsv"), emit: geno_per_cell
  path("mutations.tsv"), emit: muts
        
  
  script:
  """
  # genotype
  genotype_mutations.R \\
    --gene_bed ${gene_bed} \\
    --mutations ${mutations} \\
    --celltypes ${celltypes} \\
    --min_MQ ${params.min_MQ} \\
    --min_BQ ${params.min_BQ}
  """
}

// plot coverage
process plot_coverage {
  tag "${meta.id}_${meta.gene}"
  label 'normal10gb'
  maxRetries 10
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/", mode: "copy"
  errorStrategy 'ignore'

  input:
  tuple val(meta),
        path(gene_bed), path(gene_regions), path(gene_features),
        path(cov_per_ct), path(geno_per_ct)
  
  output:
  path("${meta.id}_${meta.gene}_*_plot.${params.plot_device}")
  path("${meta.id}_${meta.gene}_pie_mutations_plot.rds"), optional: true

  script:
  """
  plot_coverage.R \\
    --meta_id ${meta.id} \\
    --gene ${meta.gene} \\
    --gene_bed ${gene_bed} \\
    --gene_regions ${gene_regions} \\
    --gene_features ${gene_features} \\
    --cov_per_ct ${cov_per_ct} \\
    --geno_per_ct ${geno_per_ct} \\
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

  // get metadata + bam paths  
  Channel.fromPath(params.samplesheet, checkIfExists: true)
  | splitCsv(header: true)
  | map { row ->
      [[id: row.id],
       file(row.bam, checkIfExists: true),
       file(row.mutations, checkIfExists: true),
       file(row.celltypes, checkIfExists: true),
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

  // get refcds
  refcds = file(params.refcds)

  // get gene coords
  ch_gene = Channel.fromPath(params.genes).splitText().map{it -> it.trim()}
  get_gene_coords(ch_gene, gencode_gff3, refcds)

  // output bed of all gene coords
  get_gene_coords.out.gene_bed.collect()
  | collate_gene_coords

  // combine gene-x-id, subset bams
  ch_checked_bam.combine(get_gene_coords.out.gene_features)
  | map { meta, bam, bai, mutations, celltypes, cell_barcodes,
          gene, gene_bed, gene_features, exonic_bed ->
          [meta + [gene: gene],
           bam, bai, mutations, celltypes, cell_barcodes, 
           gene_bed, gene_features, exonic_bed]
  }
  | subset_bam_to_genes
  | set { ch_bam_x_gene }

  // get coverage per cell or from bulk
  
  if (params.bulk) {

    // get coverage
    bulk(ch_bam_x_gene)

  } else {

    // get coverage per celltype
    single_cell(ch_bam_x_gene)
    ch_cov_per_ct = single_cell.out.cov_per_ct
    ch_cell_bams = single_cell.out.cell_bams

    // genotype cells
    genotype_mutations(ch_cell_bams)

    // plot coverage
    plot_coverage(ch_cov_per_ct.join(genotype_mutations.out.geno_per_ct))

  }

  // // plot 5' drop-off per gene
  // ch_exonic_cov_per_gene = get_coverage_per_celltype.out.exonic_coverage.groupTuple()
  // plot_5_prime_dropoff(ch_exonic_cov_per_gene)
  
}