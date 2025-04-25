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
include { get_irods_bam  } from './modules/local/get_irods_bam.nf'
include { get_local_bam  } from './modules/local/get_local_bam.nf'
include { samtools_index } from './modules/local/samtools_index.nf'
include { MOSDEPTH       } from './modules/nf-core/mosdepth/main'
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// optionally extract the cell barcodes from the BAM
process get_cell_barcodes {
  tag "${meta.id}"
  label 'normal'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path(bam), path(bai), path("cell_barcodes.txt")

  script:
  """
  module load samtools-1.19/python-3.12.0

  # get unique cell barcodes in the BAM, trim CB:Z: prefix
  samtools view ${bam} | cut -f12- | tr "\\t" "\\n" |
  grep "CB:Z:" | sed 's/^CB:Z://g' | awk '!x[\$0]++' \
  > cell_barcodes.txt
  """
}

// optionally generate pseudo celltypes file
process fake_celltypes {
  tag "${meta.id}"
  label 'normal'

  input:
  tuple val(meta), path(bam), path(bai), path(cell_barcodes)

  output:
  tuple val(meta), path(bam), path(bai), path(cell_barcodes),
        path("celltypes.csv")

  script:
  """
  echo "barcode,celltype" > celltypes.csv
  cat cell_barcodes.txt | awk '{print \$1",all"}' >> celltypes.csv
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
        path(cell_barcodes), path(celltypes), path(mutations),
        path(gene_bed), path(gene_regions), path(gene_features)

  output:
  tuple val(meta),
        path("${meta.id}_subset_genes.bam"),
        path("${meta.id}_subset_genes.bam.bai"),
        path(cell_barcodes), path(celltypes), path(mutations),
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

// subset the BAM to only reads from barcodes of interest
process subset_bam_to_barcodes {
  tag "${meta.id}_${meta.gene}"
  label 'normal'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy",
    pattern: "*_subset.{bam,bam.bai}"
  errorStrategy 'ignore'

  input:
  tuple val(meta),
        path(bam), path(bai),
        path(cell_barcodes), path(celltypes), path(mutations),
        path(gene_bed), path(gene_regions), path(gene_features)

  output:
  tuple val(meta),
        path("${meta.id}_${meta.gene}_subset.bam"),
        path("${meta.id}_${meta.gene}_subset.bam.bai"),
        path(cell_barcodes), path(celltypes), path(mutations),
        path(gene_bed), path(gene_regions), path(gene_features)

  script:
  """
  module load samtools-1.19/python-3.12.0 

  # subset to cell barcodes
  subset-bam \
    -b ${bam} \
    --cell-barcodes ${cell_barcodes} \
    --out-bam ${meta.id}_${meta.gene}_subset.bam
  
  # index
  samtools index -@ ${task.cpus} ${meta.id}_${meta.gene}_subset.bam
  """
}

// split subsetted BAMs by cell, count coverage
process get_coverage_per_cell {
  tag "${meta.id}_${meta.gene}"
  label 'week10gb'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy",
    pattern: "*_coverage_per_cell.tsv"
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(bam), path(bai),
        path(cell_barcodes), path(celltypes), path(mutations),
        path(gene_bed), path(gene_regions), path(gene_features)

  output:
  tuple val(meta), path("${meta.id}_${meta.gene}_coverage_per_cell.tsv"),
        path(celltypes), path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features),
        emit: cov
  tuple val(meta),
        path(celltypes), path(mutations), path(gene_bed),
        path("cell_bams/*"),
        emit: cell_bams

  script:
  """
  module load samtools-1.19/python-3.12.0

  # get gene info
  chr=\$(cut -f1 ${gene_bed})
  start=\$(cut -f2 ${gene_bed})
  end=\$(cut -f3 ${gene_bed})
  echo "${meta.gene} - chr: \$chr, start: \$start, end: \$end"

  echo "getting cell barcodes"
  # get unique cell barcodes in the BAM, trim CB:Z: prefix
  samtools view ${bam} | cut -f12- | tr "\\t" "\\n" |
  grep "CB:Z:" | sed 's/^CB:Z://g' | awk '!x[\$0]++' \
  > cell_barcodes.txt

  echo "creating a bam per cell"
  rm -rf cell_bams ; mkdir -p cell_bams
  while read -r CB ; do
    # subset
    echo \$CB > \$CB.txt
    subset-bam \
      -b ${bam} \
      --cell-barcodes \$CB.txt \
      --out-bam cell_bams/\$CB.bam
    # index
    samtools index -@ ${task.cpus} cell_bams/\$CB.bam
    rm \$CB.txt
  done < cell_barcodes.txt

  # initialise gene coverage file
  rm -rf ${meta.gene} ; mkdir ${meta.gene}
  touch ${meta.gene}_cov.tsv

  echo "calculating coverage per base per cell"
  while read -r CB ; do
    (
      echo "\$CB" ;
      samtools depth \\
        -a --min-BQ ${params.min_BQ} --min-MQ ${params.min_MQ} \\
        -r \$chr:\$start-\$end cell_bams/\$CB.bam \\
      | cut -f3 ;
    ) > ${meta.gene}/\$CB.tsv
  done < cell_barcodes.txt

  echo "creating coverage file"
  (
    echo -e "chr\\tpos\\tgene" ;
    echo -e "\$chr\\t\$start\\t\$end" \
    | awk -F"\\t" -v gene=${meta.gene} -v OFS="\\t" \
      '{for(i=\$2; i<=\$3; i++) print \$1, i, gene}' ;
  ) > ${meta.gene}_coords.tsv

  echo "adding coverage per cell to coverage file"
  # (sort for consistent ordering)
  find ${meta.gene} -name "*.tsv" | sort | xargs paste \\
  > ${meta.gene}_cov.tsv

  echo "combining outputs"
  paste ${meta.gene}_coords.tsv ${meta.gene}_cov.tsv > ${meta.gene}_cov_per_cell.tsv
  ls *_cov_per_cell.tsv | head -1 | xargs head -1 \\
  > ${meta.id}_${meta.gene}_coverage_per_cell.tsv
  cat *_cov_per_cell.tsv | grep -vP "^chr\\tpos" \\
  >> ${meta.id}_${meta.gene}_coverage_per_cell.tsv
  """
}

// summarise coverage per celltype
process get_coverage_per_celltype {
  tag "${meta.id}_${meta.gene}"
  label 'long20gb'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy",
    pattern: "*_coverage_per_*.tsv"
  maxRetries 10
  
  input:
  tuple val(meta),
        path(cov),
        path(celltypes), path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features)

  output:
  tuple val(meta),
        path(gene_bed), path(gene_regions), path(gene_features),
        path("${meta.id}_${meta.gene}_coverage_per_celltype.tsv"),
        emit: cov_per_ct
  tuple val(meta),
        path("${meta.id}_${meta.gene}_coverage_per_exonic_position.tsv"),
        emit: exonic_coverage
  path("${meta.id}_${meta.gene}_coverage_per_ccds_position.tsv")

  script:
  """
  get_coverage_per_celltype.R \\
    --gene ${meta.gene} \\
    --regions ${gene_regions} \\
    --cov ${cov} \\
    --cell_barcodes ${cell_barcodes} \\
    --celltypes ${celltypes} \\
    --meta_id ${meta.id}
  """
}

// genotype the mutations
process genotype_mutations {
  tag "${meta.id}_${meta.gene}"
  label 'normal10gb'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy"
  errorStrategy 'ignore'
  
  input:
  tuple val(meta),
        path(celltypes), path(mutations),
        path(gene_bed),
        path(cell_bams, stageAs: "cell_bams/*")
  
  output:
  tuple val(meta),
        path("${meta.id}_${meta.gene}_genotyped_mutations_per_celltype.tsv"), emit: geno_per_ct
  path("${meta.id}_${meta.gene}_genotyped_mutations_per_cell.tsv"), emit: geno_per_cell
  path("mutations.tsv"), emit: muts
        
  script:
  """
  # genotype
  genotype_mutations.R \\
    --id ${meta.id} \\
    --gene ${meta.gene} \\
    --gene_bed ${gene_bed} \\
    --mutations ${mutations} \\
    --celltypes ${celltypes} \\
    --min_MQ ${params.min_MQ} \\
    --min_BQ ${params.min_BQ}
  """
}

// plot coverage
process plot_coverage_and_genotyping {
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
  tuple val(meta),
        path("${meta.id}_${meta.gene}_*_plot.${params.plot_device}"),
        emit: plots
  path("${meta.id}_${meta.gene}_pie_mutations_plot.rds"), optional: true

  script:
  """
  plot_coverage_and_genotyping.R \\
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

// knit report on mutations and coverage
process report {
  label 'normal'
  publishDir "${params.out_dir}/", mode: "copy"

  input:
  tuple val(x), val(ids), val(genes), path(plots)
  path report_rmd
      
  output:
  path("report.html")

  script:
  def c_ids = ids instanceof List ? 'c("' + ids.join('", "') + '")' : "c(\"${ids}\")"
  def c_genes = 'c("' + genes.join('", "') + '")'
  """
  #!/usr/bin/env Rscript
  rmarkdown::render(
    "${report_rmd}",
    output_file = "report.html",
    params = list(ids = ${c_ids},
                  genes = ${c_genes}))
  """
}

// plot 5 prime dropoff
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
  in_bam = \
    Channel
    .fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
    | map { meta, bam, mutations, celltypes, cell_barcodes ->
            [meta, bam]
    }

  // get input cell_barcodes
  in_cell_barcodes = \
    Channel
    .fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
    | map { meta, bam, mutations, celltypes, cell_barcodes ->
            [meta, cell_barcodes]
    }

  // get input celltypes
  in_celltypes = \
    Channel
    .fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
    | map { meta, bam, mutations, celltypes, cell_barcodes ->
            [meta, celltypes]
    }

  // get samplesheet
  in_mutations = \
    Channel
    .fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
    | map { meta, bam, mutations, celltypes, cell_barcodes ->
            [meta, mutations]
    }

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

  // get cell barcodes if not provided
  if ( params.no_cell_barcodes ) {
    ch_barcodes = get_cell_barcodes(samtools_index.out)
  } else {
    ch_barcodes = samtools_index.out.join(in_cell_barcodes)
  }

  // generate pseudo celltypes if not provided ("all")
  if ( params.no_celltypes ) {
    ch_celltypes = fake_celltypes(ch_barcodes)
  } else {
    ch_celltypes = ch_barcodes.join(in_celltypes)
  }
  
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
  ch_bam_x_gene = \
    ch_celltypes \
    .join(in_mutations) \
    .combine(get_gene_coords.out.gene_features)
    | map { meta, bam, bai, cell_barcodes, celltypes, mutations,
            gene, gene_bed, gene_features, exonic_bed ->
            [meta + [gene: gene],
              bam, bai, cell_barcodes, celltypes, mutations, 
              gene_bed, gene_features, exonic_bed]
    }
    | subset_bam_to_genes

  // get per-cell bams
  subset_bam_to_barcodes(ch_bam_x_gene)

  // get coverage per cell
  get_coverage_per_cell(subset_bam_to_barcodes.out)

  // summarise coverage
  get_coverage_per_celltype(get_coverage_per_cell.out.cov)

  // genotype mutations
  genotype_mutations(get_coverage_per_cell.out.cell_bams)

  // plot coverage and mutations
  plot_coverage_and_genotyping(get_coverage_per_celltype.out.cov_per_ct.join(genotype_mutations.out.geno_per_ct))
  
  // knit coverage plots
  report_rmd = file("${projectDir}/assets/report.Rmd", checkIfExists: true)
  ch_plots = \
    plot_coverage_and_genotyping.out.plots
    | map { meta, plots -> [meta.id, meta.gene, plots] }
    | groupTuple()
    | map { ids, genes, plots -> ['report', ids, genes, plots]}
    | map { x, ids, genes, plots_nested ->
            def plots_flat = plots_nested.flatten()
            return [x, ids, genes, plots_flat]}
  report(ch_plots, report_rmd)
  
}